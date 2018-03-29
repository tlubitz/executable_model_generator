#!/usr/bin/env python
import exportCDSBML
import sys
import copy
import molecule
import re
import SBtab
import validatorSBtab
from molecule import Molecule, MoleculeDef, Complex, Rule
from anytree import Node, RenderTree


def prepareModel(file_name: str):
    '''
    prepares the python model for the conversion to rxncon format
    '''
    ########################
    # List of subreactions #
    ########################

    def complex_structure(s_complex: Complex, structure: list=False):
        '''
        construct the simple structure of a complex in form of a list of
        elements to compare with other complexes
        '''
        global structures
        structure = ['<' + s_complex.get_name() + '>']
        for entity in s_complex.get_molecules():
            if entity.get_class().lower() == 'complex':
                ss = complex_structure(entity, structure)
                for s in ss: structure.append(s)
            else:
                try: structure.append(entity.get_protein_ref())
                except: structure.append(entity.get_name())

        return sorted(structure)

    def find_components(element):
        '''
        extracts all single components from a complex or single element
        '''
        components = []

        for el in element:
            if el.get_class().lower() == 'complex':
                comps = find_components(el.get_molecules())
                for comp in comps: components.append(comp)
            elif el.get_kmer() == '2':
                if '@' not in el.get_name():
                    k1 = copy.deepcopy(el)
                    k2 = copy.deepcopy(el)
                    k1.set_name(el.get_name() + '@' + str(el.get_index()))
                    k2.set_name(el.get_name() + '@' + str(el.get_index() + 1))
                    components.append(k1)
                    components.append(k2)
                else:
                    components.append(el)
                    components.append(el)
            else:
                components.append(el)
        return components

    def traverse(s_complex: Complex, entity: Complex, parent_node_name: Node):
        '''
        function for traversing tree structure to build the complex topology
        '''
        global rows
        global complex_ids
        global complex_structures
        global ppi_duplets

        if entity.get_id() in complex_ids: return 0
        else: complex_ids.append(entity.get_id())
        cs = '-'.join(complex_structure(entity))
        if cs in complex_structures: return 0
        else:
            entity.set_internal_id(cs)
            complex_structures.append(cs)

        complex_node = Node(entity.get_name(), parent=parent_node_name)
        leaves_local = []

        for next_leaf in entity.get_molecules():
            if next_leaf.get_class().lower() != 'complex':
                if int(next_leaf.get_kmer()) > 1:
                    create_dimer(entity, next_leaf, complex_node)
                else:
                    r_name = next_leaf.get_name()
                    leaves_local.append(next_leaf)
            else:
                equivalencies = create_equivalencies(next_leaf)
                row = '<%s>,AND,<%s>%s,,,,%s\n' % (entity.get_name().strip(),
                                                   next_leaf.get_name().strip(),
                                                   equivalencies,
                                                   entity.get_internal_id())
                if row not in rows: rows.append(row)
                traverse(entity, next_leaf, complex_node)

        if len(leaves_local) == 1:
            new_node_name = Node(leaves_local[0].get_name(),
                                 parent=complex_node)
            row = '<%s>,AND,%s@%s,,,,%s\n' % (entity.get_name().strip(),
                                              leaves_local[0].get_name(),
                                              leaves_local[0].get_index(),
                                              entity.get_internal_id())
            if row not in rows: rows.append(row)
        elif len(leaves_local) == 2:
            for i, leaf in enumerate(leaves_local):
                new_node_name = Node(leaf.get_name(), parent=complex_node)

            sortee = sorted([leaves_local[0], leaves_local[1]],
                            key=lambda x: x.get_name())
            row = '<%s>,AND,%s@%s--%s@%s,,,,%s\n' % (entity.get_name().strip(),
                                                     str(sortee[0].get_name()),
                                                     str(sortee[0].get_index()),
                                                     str(sortee[1].get_name()),
                                                     str(sortee[1].get_index()),
                                                     entity.get_internal_id())
            if row not in rows:
                rows.append(row)
                ppi_duplets.append([str(sortee[0].get_name()),
                                    str(sortee[1].get_name())])
        elif len(leaves_local) > 2:
            for i, leaf in enumerate(leaves_local):
                new_node_name = Node(leaf.get_name(), parent=complex_node)
                row = '<%s>,AND,%s@%s,,,,%s\n' % (entity.get_name().strip(),
                                                  leaf.get_name(),
                                                  leaf.get_index(),
                                                  entity.get_internal_id())
                if row not in rows: rows.append(row)

    ########################
    def create_dimer(entity: Complex, molecule: Molecule, parent_name: Node):
        '''
        function that creates a dimer
        '''
        global rows
        global complex2ppis
        global dimers
        global ppi_duplets

        dimername = molecule.get_name().strip()
        (first, second) = get_dimer_positions(entity, dimername)

        row = '<%s>,AND,<%s-Dimer>#%s@%s=%s@1#%s@%s=%s@2,,,,%s\n' % (entity.get_name().strip(),
                                                                     dimername,
                                                                     dimername,
                                                                     first,
                                                                     dimername,
                                                                     dimername,
                                                                     second,
                                                                     dimername,
                                                                     entity.get_internal_id())
        if row not in rows: rows.append(row)

        if molecule.get_name() in dimers: return 0
        else: dimers.append(molecule.get_name())
        dimer_int_id = '<%s-Dimer>-%s@1-%s@2' % (molecule.get_name().strip(),
                                                 molecule.get_name().strip(),
                                                 molecule.get_name().strip())
        row = '<%s-Dimer>,AND,%s@1--%s@2,,,,%s\n' % (molecule.get_name().strip(),
                                                     molecule.get_name().strip(),
                                                     molecule.get_name().strip(),
                                                     dimer_int_id)
        if row not in rows:
            rows.append(row)
            ppi_duplets.append([molecule.get_name().strip(),
                                molecule.get_name().strip()])
            entity.set_ppi('%s@1--%s@2' % (molecule.get_name().strip(),
                                           molecule.get_name().strip()))
            complex2ppis[parent_name.name] = ['%s@1--%s@2' % (molecule.get_name().strip(),
                                                              molecule.get_name().strip())]

    def get_dimer_positions(entity, dimername):
        '''
        retrieves the individual positions for the two dimer elements in a
        complex; this is for equivalences, e.g. #Ste12@3=Ste12@1
        '''
        first = None
        second = None
        components = find_components(entity.get_molecules())
        components = set_correct_dimer_positions(components)
        for element in components:
            try: name = element.get_name().split('@')[0]
            except: name = element.get_name()
            if first is None:
                if name == dimername:
                    first = element.get_position() + 1
                    element.set_name(name + '@' + str(element.get_position()))
            else:
                if name == dimername:
                    second = element.get_position() + 1
                    element.set_name(name + '@' + str(element.get_position()))
        return (first, second)

    def set_correct_dimer_positions(components):
        '''
        sets the correct dimer positions
        '''
        for i, c in enumerate(components): c.set_position(i)
        return components

    def fromComplexToMolecule(s_complex: Complex):
        '''
        for the rare case that we have a complex that only consists of one
        child molecule: make the complex a molecule instead
        '''
        # make new mol definition if provided
        try:
            complex_def = s_complex.get_definition()
            new_mol_def = MoleculeDef(s_complex.get_id(), s_complex.get_name())
            new_mol_def.set_notes(complex_def.get_notes())
            new_mol_def.set_species_identity(complex_def.get_species_identity())
            new_mol_def.set_protein_reference(complex_def.get_protein_reference())
            for elem in complex_def.get_modification_residues():
                new_mol_def.add_modification_residue(elem[0], elem[1])
            for elem in complex_def.get_binding_domains():
                new_mol_def.set_binding_domain(elem[0], elem[1])
        except:
            new_mol_def = MoleculeDef(s_complex.get_id(), s_complex.get_name())

        # make new molecule
        new_molecule = Molecule(s_complex.get_id(), s_complex.get_name())
        new_molecule.set_class(s_complex.get_class())
        new_molecule.set_name_cd(s_complex.get_name_cd())
        new_molecule.set_position(s_complex.get_position())
        new_molecule.set_stoich(s_complex.get_stoich())
        new_molecule.set_molecule_def(new_mol_def)
        model_object.add_molecule(new_molecule)

    def attachPPI(s_complex: Complex, ppi: str):
        '''
        attach the interacting protein objects within a complex to each other
        '''
        global obj1
        global obj2
        ppi1 = re.search('(.*)--', ppi).group(1)
        ppi2 = re.search('--(.*)', ppi).group(1)
        if '@' in ppi1: ppi1 = ppi1.split('@')[0]
        if '@' in ppi2: ppi2 = ppi2.split('@')[0]

        for entity in s_complex.get_molecules():
            if entity.get_class().lower() == 'complex':
                if len(entity.get_molecules()) == 2 and just_proteins(entity):
                    sortee = sorted([entity.get_molecules()[0],
                                     entity.get_molecules()[1]],
                                    key=lambda x: x.get_name())
                    s_complex.set_ppi('%s@%s--%s@%s' % (sortee[0].get_name(),
                                                        sortee[0].get_index(),
                                                        sortee[1].get_name(),
                                                        sortee[1].get_index()))
                else: attachPPI(entity, ppi)

    def just_proteins(s_complex: Complex):
        '''
        check if the entities of a complex are only two proteins (which implies
        that they are bonded)
        '''
        for entity in s_complex.get_molecules():
            if not entity.get_class().lower() == 'protein': return False
            if int(entity.get_kmer()) > 1: return False
        return True

    def create_equivalencies(component):
        '''
        creates the equivalencies for the complex structure contingencies,
        e.g. #Ste11@1=Ste11@4
        '''
        equivalencies = []

        for i, element in enumerate(component.get_molecules()):
            if element.get_class().lower() == 'complex':
                equivalencies.append(create_equivalencies(element))
            else:
                if not element.get_index_int(): element.set_index_int(i+1)
                try:
                    equivalencies.append('#%s@%s=%s@%s' % (element.get_name(),
                                                           element.get_index(),
                                                           element.get_name(),
                                                           element.get_index_int()))
                except:
                    print('Couldnt create an equivalency for %s in complex'\
                          '%s.' % (element.get_name(),
                                   component.get_name()))

        return ''.join(equivalencies)
    
    ########################
    # End of subreactions ##
    ########################
    # 1: get the python object model (A1)
    model_object = import_model(file_name)
    assert model_object, 'Model object could not be imported.'

    # 2: duplicate rules with more than one modifier (A2)
    modifier_count = 0
    remove_rules = []
    rules_before = len(model_object.get_rules())
    for rule in model_object.get_rules():
        if len(rule.get_modifiers()) > 1:
            for i, mod in enumerate(rule.get_modifiers()):
                new_rule = copy.copy(rule)
                new_id = rule.get_id() + '_' + str(i)
                new_rule.set_id(new_id)
                new_rule.modifiers = []
                new_rule.set_modifier(mod)
                model_object.add_rule(new_rule)
                modifier_count += 1
            remove_rules.append(rule.get_id())

    for rule in remove_rules:
        model_object.remove_rule(rule)
        rules_before -= 1

    assert len(model_object.get_rules()) == rules_before + modifier_count, 'The model rules could not be split correctly.'

    # 3: complex topology
    # This is a big one: First, we prepare an output file for the user that
    # asks for the bindings of the entities in the complex topology
    # Next, we read in the file from the user again and extract the given
    # information to store them in the python objects
    global rows
    global complex2ppis
    global complex_ids
    global complex_structures
    global dimers
    global ppi_duplets
    dimers = []
    complex2ppis = {}
    rows = []
    class2rxncon = {'gene': 'Gene', 'rna': 'mRNA'}
    complex_ids = []
    ppi_duplets = []
    complex_structures = []

    # All functions assembled, now starting to build tree and output file
    root = Node(file_name[:-4])
    for s_complex in model_object.get_complexes():
        if s_complex.get_id() in complex_ids: continue
        else: complex_ids.append(s_complex.get_id())

        cs = '-'.join(complex_structure(s_complex))
        if cs in complex_structures: continue
        else:
            s_complex.set_internal_id(cs)
            complex_structures.append(cs)

        # collect the leaves
        complex_node = Node(s_complex.get_name(), parent=root)
        leaves = []
        for entity in s_complex.get_molecules():
            # for the rare case that we have a complex that only consists
            # of one child molecule: make the complex a molecule instead
            if len(s_complex.get_molecules()) == 1 and entity.get_class().lower() != 'complex':
                fromComplexToMolecule(s_complex)
                try: model_object.remove_complex(s_complex.get_id())
                except: pass
                break
            if entity.get_class().lower() != 'complex':
                if int(entity.get_kmer()) > 1:
                    create_dimer(s_complex, entity, complex_node)
                else: leaves.append(entity)
            else:
                equivalencies = create_equivalencies(entity)
                row = '<%s>,AND,<%s>%s,,,,%s\n' % (s_complex.get_name().strip(),
                                                   entity.get_name().strip(),
                                                   equivalencies,
                                                   s_complex.get_internal_id())
                if row not in rows: rows.append(row)
                traverse(s_complex, entity, complex_node)

        # write leaf content to anytree (for proof check mainly) and to output
        # rows. If only 1 leaf:
        if len(leaves) == 1:
            r_name = leaves[0].get_name()
            leaf_node_name = Node(r_name, parent=complex_node)
            row = '<%s>,AND,%s@%s,,,,%s\n' % (s_complex.get_name().strip(),
                                              r_name, leaves[0].get_index(),
                                              s_complex.get_internal_id())
            if row not in rows: rows.append(row)
        # if 2 leaves
        elif len(leaves) == 2:
            # if complex consists of more than 2 leaves (e.g. also complexes)
            if len(s_complex.get_molecules()) != 2:
                for i, leaf in enumerate(leaves):
                    r_name = leaf.get_name().strip()
                    leaf_node_name = Node(r_name, parent=complex_node)
                    row = '<%s>,AND,%s@%s,,,,%s\n' % (s_complex.get_name().strip(),
                                                      r_name.strip(),
                                                      leaf.get_index(),
                                                      s_complex.get_internal_id())
                    if row not in rows: rows.append(row)
            # if complex ONLY consists of 2 leaves
            else:
                sortee = sorted([leaves[0], leaves[1]],
                                key=lambda x: x.get_name())
                row = '<%s>,AND,%s@%s--%s@%s,,,,%s\n' % (s_complex.get_name().strip(),
                                                         sortee[0].get_name(),
                                                         sortee[0].get_index(),
                                                         sortee[1].get_name(),
                                                         sortee[1].get_index(),
                                                         s_complex.get_internal_id())
                if row not in rows:
                    rows.append(row)
                    ppi_duplets.append([sortee[0].get_name(),
                                        sortee[1].get_name()])
                    s_complex.set_ppi('%s@%s--%s@%s' % (sortee[0].get_name(),
                                                        sortee[0].get_index(),
                                                        sortee[1].get_name(),
                                                        sortee[1].get_index()))
                for leaf in leaves:
                    leaf_node_name = Node(leaf.get_name(), parent=complex_node)
        # if more than 2 leaves
        elif len(leaves) > 2:
            for i, leaf in enumerate(leaves):
                r_name = leaf.get_name()
                leaf_node_name = Node(r_name, parent=complex_node)
                row = '<%s>,AND,%s@%s,,,,%s\n' % (s_complex.get_name().strip(),
                                                  r_name.strip(),
                                                  leaf.get_index(),
                                                  s_complex.get_internal_id())
                if row not in rows: rows.append(row)

    # Output of a tree to commandline and a csv file for the user
    # from anytree.dotexport import RenderTreeGraph
    # RenderTreeGraph(root_objects).to_picture("tree.png")
    for pre, fill, node in RenderTree(root):
        # print("%s%s" % (pre, node.name))
        pass

    sbtab_string = '!!SBtab TableType="rxnconContingencyList" '\
                   'TableName="ContingencyList" Document="%s" '\
                   '\n' % (file_name[:-4])
    sbtab_string += '!UID:Contingency,!Target,!Contingency,!Modifier,'\
                    '!Reference:Identifiers:pubmed,!Quality,!Comment,'\
                    '!InternalComplexID\n'
    for i, row in enumerate(sorted(rows)):
        sbtab_string += str(i + 1) + ',' + row


    ff = re.search('(.*)/(.*)', file_name[:-4].lower()).group(2)
    filename = 'files/%s_bindings.csv' % (ff)
    sbtab_contingency = SBtab.SBtabTable(sbtab_string, filename)
    sbtab_contingency_valid = validatorSBtab.ValidateTable(sbtab_contingency,
                                                           filename)
    warnings = sbtab_contingency_valid.return_output()
    if warnings != []:
        print('Warnings for file %s:\n' % filename)
        for warning in warnings:
            print(warning)
    sbtab_contingency.write(filename)         
    

    # 3.2 Now reading the file back in and storing the provided information
    print('The file bindings_filename.csv has been written to directory inputoutput.'\
          'Please open the file, add the requested bindings of the complex'\
          'topology, and save it. If you"re done, please hit enter and this'\
          'script will proceed.')
    input('')
    get_file = False

    while get_file is False:
        try:
            filename_re = filename[:-4] + '_re.csv'
            domain_file = open(filename_re, 'r').read()
            get_file = True
        except:
            print('The domain file could not be read. Please make sure the'\
                  'file name and location have not changed and then hit en'\
                  'ter again.')
            input('')

    sbtab_contingency_re = SBtab.SBtabTable(domain_file, filename_re)
    sbtab_contingency_re_valid = validatorSBtab.ValidateTable(sbtab_contingency_re,
                                                              filename_re)
    warnings = sbtab_contingency_re_valid.return_output()
    if warnings != []:
        print('Warnings for file %s:\n' % filename_re)
        for warning in warnings:
            print(warning)

    saver = []
    for row in sbtab_contingency_re.value_rows:
        saver.append(','.join(row[1:]))
        row_id = row[0]
        row_c = re.search('<(.*)>', row[1]).group(1)
        if not row[3].startswith('<') and not row[3].endswith('>'):
            row_original_parts = row[3].split('--')
            row_parts_sorted = sorted(row_original_parts)
            row_ppi = '--'.join(row_parts_sorted)
        else: continue
        if row_c in complex2ppis.keys():
            ppis = complex2ppis[row_c]
            if row_ppi not in ppis: ppis.append(row_ppi)
            complex2ppis[row_c] = ppis
        else:
            complex2ppis[row_c] = [row_ppi]

    def collect_ppis(s_complex):
        '''
        collect all ppis that are attached to complex and its subcomplexes
        '''
        ppis = complex2ppis[s_complex.get_name()]

        for elem in s_complex.get_molecules():
            if elem.get_class().lower() == 'complex':
                more_ppis = collect_ppis(elem)
                for ppi in more_ppis: ppis.append(ppi)
            elif elem.get_kmer() == 2:
                try: ppis.append(complex2ppis[elem.get_name()])
                except: pass
            else:
                try: ppis.append(complex2ppis[elem])
                except: pass

        return ppis

    # Information is read from the file, now incorporate into Python objects
    for s_complex in model_object.get_complexes():
        if s_complex.get_name() in complex2ppis.keys():
            ppis = collect_ppis(s_complex)
            for ppi in ppis:
                ppi1 = re.search('(.*)--', ppi).group(1)
                ppi2 = re.search('--(.*)', ppi).group(1)
                sortee = sorted([ppi1, ppi2])
                s_complex.set_ppi('%s--%s' % (sortee[0], sortee[1]))

                for st_complex in model_object.get_complexes():
                    attachPPI(st_complex, ppi)

        # 2 proteins in 1 complex are automatically interacting
        if len(s_complex.get_molecules()) == 2 and just_proteins(s_complex):
            sortee = sorted([s_complex.get_molecules()[0],
                             s_complex.get_molecules()[1]],
                            key=lambda x: x.get_name())
            s_complex.set_ppi('%s@%s--%s@%s' % (sortee[0].get_name(),
                                                sortee[0].get_index(),
                                                sortee[1].get_name(),
                                                sortee[1].get_index()))

    for mol in model_object.get_molecules():
        if mol.get_kmer() == '2':
            try:
                ppis = complex2ppis[mol.get_name() + '-Dimer']
                for ppi in ppis:
                    ppi1 = re.search('(.*)--', ppi).group(1)
                    ppi2 = re.search('--(.*)', ppi).group(1)
                    if '@' in ppi1: ppi1 = ppi1.split('@')[0]
                    if '@' in ppi2: ppi2 = ppi2.split('@')[0]
                    sortee = sorted(ppi1, ppi2)
                    mol.set_ppi('%s--%s' % (sortee[0], sortee[1]))
                    for st_complex in model_object.get_complexes():
                        attachPPI(st_complex, ppi)
            except: pass

    # 4: create list of component states (the molecule definitions in the
    #    Python model object) (A6)
    component_states = []
    for molecule in model_object.get_molecules():
        mol_def = molecule.get_definition()
        component_states.append(mol_def)

    return model_object, ppi_duplets, saver


def import_model(file_name: str):
    '''
    imports the python model that has been generated from CD SBML
    '''
    model_object = exportCDSBML.export_CDSBML(file_name)
    return model_object


if __name__ == '__main__':

    try: sys.argv[1]
    except: sys.exit()

    file_name = sys.argv[1]
    intermediate, duplets, saver = prepareModel(file_name)
