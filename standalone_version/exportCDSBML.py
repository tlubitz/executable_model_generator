#!/usr/bin/env python
import re
import libsbml
import numpy
import sys
import copy
from molecule import Molecule, MoleculeDef, Complex, Rule
from model import Model, Parameter, InitialCondition


def export_CDSBML(sbml_file_name: str):
    '''
    Takes an SBML file with CD annotations and exports the inherent information
    into a Model object. This is done in 5 subsequent steps:
    1: parameters (if given in SBML)
    2: molecule definitions
    3: molecules (first cd then sbml species) + init conditions
    4: complexes + init conditions
    5: reactions/rules
    '''
    # function used for recursive traversion of the hierarchy tree
    def backtrace_complexes(species: libsbml.XMLNode):
        species_ident = species.getAttrValue(0)
        species_name_cd = species.getAttrValue(1)

        species_annotation = species.getChild('annotation')
        complex_species = species_annotation.getChild('complexSpecies').getChild(0).getCharacters()
        species_class = species_annotation.getChild('speciesIdentity').getChild('class').getChild(0).getCharacters()
        protein_reference = species_annotation.getChild('speciesIdentity').getChild('proteinReference').getChild(0).getCharacters()
        kmer_association = species_annotation.getChild('speciesIdentity').getChild('state').getChild('homodimer').getChild(0).getCharacters()

        try: molecules = [mol_id2molecule[species_ident]]
        except: molecules = [complex_id2complex[species_ident]]

        return molecules

    # make an instance of the SBML model object
    reader = libsbml.SBMLReader()
    sbml = reader.readSBML(sbml_file_name)
    model = sbml.getModel()

    # check if the model holds CD annotations
    if not determineCD(model): sys.exit()

    # initialise the model
    new_model = Model(model.getId(), model.getName())

    ##################################################
    # 1: transfer parameters from SBML to target object
    for parameter in model.getListOfParameters():
        if parameter.getName() != '':
            new_parameter = Parameter(parameter.getName(),
                                      parameter.getValue())
        else: new_parameter = Parameter(parameter.getId(),
                                        parameter.getValue())
        new_model.add_parameter(new_parameter)

    ##################################################
    # 2: extract CD SBML molecules & definitions
    # get CD lists
    model_extension = model.getAnnotation().getChild(0)
    list_of_cd_species = model_extension.getChild('listOfIncludedSpecies')
    list_of_proteins = model_extension.getChild('listOfProteins')

    # collect the CD proteins which correspond to molecule definitions
    definition_id2object = {}
    complex_id2complex = {}

    for i in range(list_of_proteins.getNumChildren()):
        current_protein = list_of_proteins.getChild(i)
        new_molecule_def = MoleculeDef(current_protein.getAttrValue(0),
                                       current_protein.getAttrValue(1))
        mod_residues = current_protein.getChild('listOfModificationResidues')
        for k in range(mod_residues.getNumChildren()):
            current_modification = mod_residues.getChild(k)
            res_id = current_modification.getAttrValue(1)
            try: res_name = current_modification.getAttrValue(2)
            except: continue
            res_angle = current_modification.getAttrValue(0)
            res_side = current_modification.getAttrValue(3)
            res_state = None
            if res_name == 'none' or res_name == '?': res_name = None
            if res_side == 'none' or res_side == '?': res_side = None
            if res_angle == 'none' or res_angle == '?': res_angle = None
            if res_id not in new_molecule_def.modification_residues.keys():
                new_molecule_def.add_modification_residue(res_id, res_state,
                                                          res_name, res_angle,
                                                          res_side)
        definition_id2object[current_protein.getAttrValue(0)] = new_molecule_def

    ##################################################
    # 3a: get CD species
    # this also preprocesses the objects to traverse the hierarchy tree later
    # on during complex buildup
    complex_species2object = {}
    complex_species2ident = {}
    mol_id2molecule = {}
    nucleotides = ['GDP', 'GTP', 'mRNA']
    for i in range(list_of_cd_species.getNumChildren()):
        species = list_of_cd_species.getChild(i)
        species_ident = species.getAttrValue(0)
        species_name_cd = species.getAttrValue(1).title()
        species_name_cd = re.sub(r'\W+', '', species_name_cd)
        for nucleotide in nucleotides:
            if nucleotide.title() in species_name_cd:
                species_name_cd = species_name_cd.replace(nucleotide.title(),
                                                          nucleotide)

        species_annotation = species.getChild('annotation')
        complex_species = species_annotation.getChild('complexSpecies').getChild(0).getCharacters()
        species_class = species_annotation.getChild('speciesIdentity').getChild('class').getChild(0).getCharacters()
        protein_reference = species_annotation.getChild('speciesIdentity').getChild('proteinReference').getChild(0).getCharacters()
        modifications = species_annotation.getChild('speciesIdentity').getChild('state').getChild('listOfModifications')
        kmer_association = species_annotation.getChild('speciesIdentity').getChild('state').getChild('homodimer').getChild(0).getCharacters()

        if species_class.lower() == 'complex':
            new_complex = Complex(species_ident, species_name_cd)
            new_complex.set_class(species_class)
            complex_id2complex[species_ident] = new_complex
        else:
            new_molecule = Molecule(species_ident, species_name_cd)
            if species_class != '': new_molecule.set_class(species_class)
            if species_name_cd != '': new_molecule.set_name_cd(species_name_cd)
            if protein_reference != '':
                new_molecule.set_protein_ref(protein_reference)
            if kmer_association != '': new_molecule.set_kmer(kmer_association)
            if modifications != '':
                try:
                    current_moldef = copy.deepcopy(definition_id2object[protein_reference])
                    for j in range(modifications.getNumChildren()):
                        if modifications.getChild(j).getAttrValue(0) not in current_moldef.modification_residues.keys():
                            current_moldef.add_modification_residue(residue=modifications.getChild(j).getAttrValue(0), state=modifications.getChild(j).getAttrValue(1))
                        elif modifications.getChild(j).getAttrValue(1) is not None:
                            current_moldef.set_mod_state(modifications.getChild(j).getAttrValue(0), modifications.getChild(j).getAttrValue(1))
                        else: pass
                    new_molecule.set_molecule_def(current_moldef)
                except:
                    mods = []
                    for j in range(modifications.getNumChildren()):
                        mod = (modifications.getChild(j).getAttrValue(0),
                               modifications.getChild(j).getAttrValue(1))
                        mods.append(mod)

            # if there is no proper definition (e.g. for DEGRADED),
            # then set an empty rudimentary molecule definition
            if not new_molecule.molecule_def:
                new_molecule_def = MoleculeDef(species_ident, species_name_cd)
                new_molecule.set_molecule_def(new_molecule_def)

            new_model.add_molecule(new_molecule)
            if species_class.lower() == 'simple_molecule':
                new_model.add_simple_molecule(species_name_cd)
            mol_id2molecule[species_ident] = new_molecule

        if complex_species in complex_species2object.keys():
            speciess = complex_species2object[complex_species]
            pspecies = complex_species2ident[complex_species]
            speciess.append(species)
            pspecies.append(species_ident)
            complex_species2object[complex_species] = speciess
            complex_species2ident[complex_species] = pspecies
        else:
            complex_species2object[complex_species] = [species]
            complex_species2ident[complex_species] = [species_ident]

    for entry in complex_id2complex:
        current_complex = complex_id2complex[entry]
        molecules = complex_species2ident[entry]
        for element in molecules:
            try: current_element = mol_id2molecule[element]
            except: current_element = complex_id2complex[element]
            current_complex.add_molecule(current_element)
        complex_id2complex[entry] = current_complex

    ##################################################
    # 3b: get SBML species list
    list_of_sbml_species = model.getListOfSpecies()

    for species in list_of_sbml_species:
        species_ident = species.getId()
        species_name_sbml = species.getName().title()
        species_name_sbml = re.sub(r'\W+', '', species_name_sbml)
        species_compart = species.getCompartment()
        for nucleotide in nucleotides:
            if nucleotide.title() in species_name_sbml:
                species_name_sbml = species_name_sbml.replace(nucleotide.title(),
                                                              nucleotide)

        species_extension = species.getAnnotation().getChild(0)
        species_position = species_extension.getChild('positionToCompartment').getChild(0).getCharacters()
        species_identity = species_extension.getChild('speciesIdentity')
        species_class = species_identity.getChild('class').getChild(0).getCharacters()
        species_name_cd = species_identity.getChild('name').getChild(0).getCharacters()
        kmer_association = species_identity.getChild('state').getChild('homodimer').getChild(0).getCharacters()
        protein_reference = species_identity.getChild('proteinReference').getChild(0).getCharacters()

        if species_class.lower() == 'complex': continue

        new_molecule = Molecule(species_ident, species_name_sbml)

        if species_compart != '': new_molecule.set_compartment(species_compart)
        if species_class != '':
            new_molecule.set_class(species_class)
            if species_class.lower() == 'gene':
                new_molecule.set_name(species_name_sbml + 'Gene')
            elif species_class.lower() == 'rna':
                new_molecule.set_name(species_name_sbml + 'mRNA')
        if species_name_cd != '': new_molecule.set_name_cd(species_name_cd)
        if kmer_association != '': new_molecule.set_kmer(kmer_association)

        list_of_catalyzed_reactions = species_extension.getChild('listOfCatalyzedReactions')
        if list_of_catalyzed_reactions.getNumChildren() > 0:
            catalyzed_reactions = []
            for i in range(list_of_catalyzed_reactions.getNumChildren()):
                cat_reaction_id = list_of_catalyzed_reactions.getChild(i)
            catalyzed_reactions.append(cat_reaction_id)
            new_molecule.set_cat_reactions(catalyzed_reactions)

        list_of_modifications = species_identity.getChild('state').getChild('listOfModifications')
        if list_of_modifications.getNumChildren() > 0:
            try:
                current_moldef = copy.deepcopy(definition_id2object[protein_reference])
                for i in range(list_of_modifications.getNumChildren()):
                    modification = list_of_modifications.getChild(i)
                    if modification.getAttrValue('residue') not in current_moldef.modification_residues.keys():
                        current_moldef.add_modification_residue(residue=modification.getAttrValue('residue'),
                                                                state=modification.getAttrValue('state'))
                    elif modification.getAttrValue('state') is not None:
                        current_moldef.set_mod_state(modification.getAttrValue(0),
                                                     modification.getAttrValue(1))
                    else: pass
                new_molecule.set_molecule_def(current_moldef)
            except: pass

        if not new_molecule.molecule_def:
            try:
                prot_mol = definition_id2object[protein_reference]
                new_molecule.set_molecule_def(prot_mol)
            except:
                # if there is no proper definition (e.g. for DEGRADED),
                # then set an empty rudimentary molecule definition
                new_molecule_def = MoleculeDef(species_ident,
                                               species_name_sbml)
                new_molecule.set_molecule_def(new_molecule_def)

        new_model.add_molecule(new_molecule)
        if species_class.lower() == 'simple_molecule':
            new_model.add_simple_molecule(species_name_sbml)
        mol_id2molecule[species_ident] = new_molecule

        if species.isSetInitialAmount():
            initial_value = species.getInitialAmount()
        elif species.isSetInitialConcentration():
            initial_value = species.getInitialConcentration()
        else: initial_value = False

        initial_condition = InitialCondition(new_molecule, initial_value)
        new_model.add_initial_condition(initial_condition)

    ##################################################
    # 4: now that we have parameters, molecule definitions, molecules,
    # and initial conditions, we focus on complexes
    for species in list_of_sbml_species:
        species_ident = species.getId()
        species_name_sbml = species.getName().title()
        species_name_sbml = re.sub(r'\W+', '', species_name_sbml)
        species_compart = species.getCompartment()
        for nucleotide in nucleotides:
            if nucleotide.title() in species_name_sbml:
                species_name_sbml = species_name_sbml.replace(nucleotide.title(),
                                                              nucleotide)

        species_extension = species.getAnnotation().getChild(0)
        species_position = species_extension.getChild('positionToCompartment').getChild(0).getCharacters()
        species_identity = species_extension.getChild('speciesIdentity')
        species_class = species_identity.getChild('class').getChild(0).getCharacters()
        species_name_cd = species_identity.getChild('name').getChild(0).getCharacters()
        protein_reference = species_identity.getChild('proteinReference').getChild(0).getCharacters()

        if species_class.lower() != 'complex': continue

        new_complex = Complex(species_ident, species_name_sbml)
        for involved_molecule in complex_species2object[species_ident]:
            complex_molecule = backtrace_complexes(involved_molecule)
            if complex_molecule is None: continue
            for molecule in complex_molecule:
                new_complex.add_molecule(molecule)

        if species_compart != '': new_complex.set_compartment(species_compart)
        if species_class != '': new_complex.set_class(species_class)
        if species_name_cd != '': new_complex.set_name_cd(species_name_cd)
        if protein_reference != '':
            new_complex.set_protein_ref(protein_reference)

        list_of_catalyzed_reactions = species_extension.getChild('listOfCatalyzedReactions')
        if list_of_catalyzed_reactions.getNumChildren() > 0:
            catalyzed_reactions = []
            for i in range(list_of_catalyzed_reactions.getNumChildren()):
                cat_reaction_id = list_of_catalyzed_reactions.getChild(i)
            catalyzed_reactions.append(cat_reaction_id)
            new_complex.set_cat_reactions(catalyzed_reactions)

        complex_id2complex[species_ident] = new_complex

        if species.isSetInitialAmount():
            initial_value = species.getInitialAmount()
        elif species.isSetInitialConcentration():
            initial_value = species.getInitialConcentration()
        else: initial_value = False

        initial_condition = InitialCondition(new_complex, initial_value)
        new_model.add_initial_condition(initial_condition)
        new_model.add_complex(new_complex)

    for c in new_model.get_complexes():
        c.sort_molecules()

    ###########################################################
    # 4b: the complex entities need proper indices
    def traverse(elem, c_index):
        '''
        traverses the tree of a complex (also iteratively if required)
        '''
        index_internal = 1
        if elem.get_class().lower() == 'complex': elem.sort_molecules()
        for el in elem.get_molecules():
            if el.get_class().lower() == 'complex':
                (el, c_index) = traverse(el, c_index)
            else:
                el.set_index_int(index_internal)
                el.set_index(c_index)
                if el.get_kmer() != 1:
                    index_internal += 2
                    c_index += 2
                else:
                    index_internal += 1
                    c_index += 1

        return (elem, c_index)

    for c in new_model.get_complexes():
        c.sort_molecules()
        complex_index = 1
        for elem in c.get_molecules():
            if elem.get_class().lower() == 'complex':
                (elem, complex_index) = traverse(elem, complex_index)
            else:
                elem.set_index(complex_index)
                if elem.get_kmer() != 1: complex_index += 2
                else: complex_index += 1

    ###########################################################
    # 5: translate reactions to rules of the model
    list_of_reactions = model.getListOfReactions()
    for reaction in list_of_reactions:
        reaction_ident = reaction.getId()
        reaction_name = reaction.getName()

        new_rule = Rule(reaction_ident, reaction_name)
        new_rule.set_direction(reaction.getReversible())

        extension = reaction.getAnnotation().getChild(0)
        reaction_type = extension.getChild('reactionType')
        new_rule.set_reaction_type(reaction_type.getChild(0).getCharacters())

        list_of_modification = extension.getChild('listOfModification')
        for i in range(list_of_modification.getNumChildren()):
            try: modification = list_of_modification.getChild(i).getAttrValue('type').lower()
            except: modification = False
            if not modification: continue
            try: modifier = list_of_modification.getChild(i).getAttrValue('modifiers').lower()
            except: modifier = False
            if not modifier: continue
            new_rule.set_modification_type(modifier, modification)

        for reactant in reaction.getListOfReactants():
            try: react = mol_id2molecule[reactant.getSpecies()]
            except:
                try: react = complex_id2complex[reactant.getSpecies()]
                except: react = False
            react.set_stoich(reactant.getStoichiometry())
            new_rule.set_left_hand_molecule(react)

        for product in reaction.getListOfProducts():
            try: prod = mol_id2molecule[product.getSpecies()]
            except:
                try: prod = complex_id2complex[product.getSpecies()]
                except: prod = False
            prod.set_stoich(product.getStoichiometry())
            new_rule.set_right_hand_molecule(prod)

        for modifier in reaction.getListOfModifiers():
            try: mod = mol_id2molecule[modifier.getSpecies()]
            except:
                try: mod = complex_id2complex[modifier.getSpecies()]
                except: mod = False

            new_rule.set_modifier(mod)
        new_model.add_rule(new_rule)

    return new_model


def determineCD(model: libsbml.Model):
    '''
    if the provided SBML is exported from cell designer,
    loads of meta information need to be processed.
    '''
    cd = False
    for i in range(model.getNamespaces().getNumNamespaces()):
        if model.getNamespaces().getPrefix(i) == 'celldesigner':
            cd = True

    return cd


if __name__ == '__main__':

    try: sys.argv[1]
    except: sys.exit()

    file_name = sys.argv[1]
    model_object = export_CDSBML(file_name)
    print("")
