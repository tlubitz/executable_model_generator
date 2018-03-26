#!/usr/bin/env python
import prepareConversion
import sys
import copy
import molecule
import re
from molecule import Molecule, MoleculeDef, Complex, Rule

def convert(file_name: str):
    '''
    prepares the python model for the conversion to rxncon format
    '''
    ########################
    # List of subreactions #
    ########################
    def pheno_check(rule, lhs, rhs):
        '''
        checks if the class of the products contains a phenotype, which
        needs to be handled differently than other reactions
        '''
        global events

        mods = rule.get_modification_types()
        if 'boolean_logic_gate_and' in mods.values(): bool_term = 'AND'
        else: bool_term = 'OR'

        for elem in rhs:
            if elem.get_class().lower() == 'phenotype':
                e_list = []
                if elem.get_name() in events.keys():
                    c_events = events[elem.get_name()]
                    for e in c_events: e_list.append(e)
                    if lhs[0].get_class().lower() == 'complex':
                        e_list.append('<%s_Influence>,%s,<%s>' % (elem.get_name(),
                                                                  bool_term,
                                                                  lhs[0].get_name()))
                    else:
                        e_list.append('<%s_Influence>,%s,%s' % (elem.get_name(),
                                                                bool_term,
                                                                lhs[0].get_name()))
                    events[elem.get_name()] = e_list
                else:
                    if lhs[0].get_class().lower() == 'complex':
                        contingencies = ['[%s],!,<%s_Influence>' % (elem.get_name(),
                                                                    elem.get_name()),
                                         '<%s_Influence>,%s,<%s>' % (elem.get_name(),
                                                                     bool_term,
                                                                     lhs[0].get_name())]
                    else: contingencies = ['[%s],!,<%s_Influence>' % (elem.get_name(),
                                                                      elem.get_name()),
                                           '<%s_Influence>,%s,%s' % (elem.get_name(),
                                                                     bool_term,
                                                                     lhs[0].get_name())]
                    events[elem.get_name()] = contingencies
                return True
            else: return False
        return False

    def trans_check(rule):
        '''
        checks if the rule is a translation, transcription, or similar
        '''
        global contingencies
        global r2terms
        global catalyst_replacements
        lhs = rule.get_left_hand_side()
        rhs = rule.get_right_hand_side()

        try:
            m = rule.get_modifiers()[0]
            m_id = m.get_id()
            if m.get_class().lower() == 'complex':
                components = find_components([m])
                cs = '['
                for c in components:
                    if '@' in c.get_name(): name = c.get_name().split('@')[0]
                    else: name = c.get_name()
                    cs += name + '|'
                modifier = cs[:-1] + ']'
            else: modifier = m.get_name()
        except:
            modifier = 'UnknownComponent'
            try: m_id
            except: m_id = False

        mt = rule.get_modification_types()

        if m_id in mt.keys(): modification_type = mt[m_id]
        else: modification_type = False

        if len(lhs) == 1 and len(rhs) == 1:
            lc = lhs[0].get_class().lower()
            rc = rhs[0].get_class().lower()
            if (lc, rc) == ('gene', 'rna'):
                if modification_type == 'catalysis':
                    if '_' in rule.get_id():
                        if '@' in modifier: modifier = modifier.split('@')[0]
                        if '@' in rhs[0].get_name():
                            k = rhs[0].get_name().split('@')[0]
                        else: b = rhs[0].get_name()
                        try: k = b.split('mRNA')[0] + 'Gene'
                        except: k = 'Gene'
                        term = '%s_TRSC_%s,%s,,,TRSC,%s,,,,,' % (modifier, k,
                                                                 modifier, k)
                        if '|' in modifier:
                            if term in catalyst_replacements.keys():
                                reactions = catalyst_replacements[term]
                                reactions.append(rule.get_id())
                                catalyst_replacements[term] = reactions
                            else: catalyst_replacements[term] = [rule.get_id()]
                        s_r = rule.get_id().split('_')[0]
                        if s_r in r2terms.keys():
                            k = r2terms[s_r]
                            k.append(term)
                            r2terms[s_r] = k
                        else: r2terms[s_r] = [term]
                        return term
                elif modification_type == 'inhibition':

                    try:
                        clean_id = rule.get_id().split('_')[0]
                        terms = r2terms[clean_id]
                        if len(terms) > 1:
                            term = []
                            for elem in terms:
                                t = '%s_%s,x,%s' % (elem.split(',')[0],
                                                    rule.get_id(), modifier)
                                contingencies.append(t)
                        else:
                            t = '%s_%s,x,%s' % (terms[0].split(',')[0],
                                                rule.get_id(), modifier)
                            contingencies.append(t)
                    except: return None
            elif (lc, rc) == ('rna', 'protein'):
                if '@' in modifier: modifier = modifier.split('@')[0]
                if '@' in rhs[0].get_name():
                    k = rhs[0].get_name().split('@')[0]
                else: k = rhs[0].get_name()
                term = '%s_TRSL_%s,%s,,,TRSL,%s,,,,,' % (modifier, k,
                                                         modifier, k)
                if '|' in modifier:
                    if term in catalyst_replacements.keys():
                        reactions = catalyst_replacements[term]
                        reactions.append(rule.get_id())
                        catalyst_replacements[term] = reactions
                    else: catalyst_replacements[term] = [rule.get_id()]
                return term
            # This is transport: but we omit this for now
            # elif ln == rn and l_loc is not None and r_loc is not None:
            #     if l_loc != r_loc:
            #         term = '%s,TRANSPORT,%s' % (modifier, ln)
            #         return term
            #
            else: return None

    def find_components(element):
        '''
        extracts all single components from a complex or single element
        '''
        components = []
        last_index = 1

        for el in element:
            if el.get_class().lower() == 'complex':
                comps = find_components(el.get_molecules())
                for comp in comps: components.append(comp)
                last_index = comp.get_index()
            elif el.get_kmer() == '2':
                if el.get_index() is None: el.set_index(last_index)
                k1 = copy.deepcopy(el)
                k2 = copy.deepcopy(el)
                k2.set_index(k1.get_index() + 1)
                components.append(k1)
                components.append(k2)
                last_index = el.get_index()
            else:
                if el.get_index() is None: el.set_index(last_index)
                components.append(el)

        return components

    def delta_c(rule):
        '''
        evaluates changes in component composition of left and right hand side
        of the rule
        '''
        r_new = []
        l_deg = []
        r = [p.get_name() for p in rule.get_single_products()]
        l = [p.get_name() for p in rule.get_single_substrates()]

        for elem in rule.get_single_substrates():
            if elem.get_name() not in r:
                l_deg.append(elem)
                rule.set_deg(elem)

        for elem in rule.get_single_products():
            if elem.get_name() not in l:
                r_new.append(elem)
                rule.set_synth(elem)

        if r_new != []: rule.rxn_type = 'synthesis'
        if l_deg != []: rule.rxn_type = 'degradation'
        if r_new != [] and l_deg != []: rule.rxn_type = 'synthesis/degradation'

    def delta_ppi(lhs_ppis, rhs_ppis):
        '''
        evaluates changes in ppi composition of left and right hand side
        '''
        global log_file
        log_file += 'PPIs LHS: %s\n' % (lhs_ppis)
        log_file += 'PPIs RHS: %s\n' % (rhs_ppis)
        lost = []
        new = []

        if sorted(lhs_ppis) != sorted(rhs_ppis):
            lost = set(lhs_ppis) - set(rhs_ppis)
            new = set(rhs_ppis) - set(lhs_ppis)

        return [lost, new]

    def delta_m(rule, lhs, rhs):
        '''
        evaluates changes in modification states of left and right hand side
        '''
        global log_file
        global writelines
        global catalyst_replacements
        global mod_triplets
        global oppress_empty_residues

        pls = {'phosphorylated': 'P+', 'acetylated': 'Ace+',
               'ubiquitinated': 'Ub+', 'palmytoylated': 'palm+',
               'prenylated': 'pre+', 'autophosphorylated': 'AP+'}
        mns = {'phosphorylated': 'P-', 'acetylated': 'Ace-',
               'ubiquitinated': 'Ub-', 'palmytoylated': 'palm-',
               'prenylated': 'pre-', 'autophosphorylated': 'AP-'}
        try:
            modifier = rule.get_modifiers()[0]
            if modifier.get_class().lower() == 'complex':
                components = find_components([modifier])
                cs = '['
                for c in components:
                    if '@' in c.get_name(): name = c.get_name().split('@')[0]
                    else: name = c.get_name()
                    cs += name + '|'
                modifier = cs[:-1] + ']'
            else:
                modifier = modifier.get_name()
                if '@' in modifier:
                    modifier = modifier.get_name().split('@')[0]
        except: modifier = 'UnknownDemodifier'

        pairs = find_pairs(lhs, rhs)
        terms = []

        for pair in pairs:
            mods_l = get_mods(pair[0])
            mods_r = get_mods(pair[1])
            log_file += 'Modifications LHS: %s\n' % (mods_l)
            log_file += 'Modifications RHS: %s\n' % (mods_r)
            for component in mods_l:
                # try:
                    try:
                        residues_l = mods_l[component]
                        residues_r = mods_r[component]
                    except: continue
                    for residue in residues_l:
                        res_l = residues_l[residue]
                        res_r = residues_r[residue]
                        if res_l[1] is None and res_r[1] is not None:
                            try: rxn = pls[res_r[1].lower()]
                            except: rxn = 'Mod+'
                            rule.mod_type = rxn
                            res = res_r[0]
                            if res is None: res = residue
                            if '@' in component:
                                component = component.split('@')[0]
                            termx = '%s_%s_%s_[(%s)],%s,,,%s,%s,,%s,,,' % (modifier,
                                                                           rxn,
                                                                           component,
                                                                           res,
                                                                           modifier,
                                                                           rxn,
                                                                           component,
                                                                           res)
                            if '|' in modifier:
                                if termx in catalyst_replacements.keys():
                                    reactions = catalyst_replacements[termx]
                                    reactions.append(rule.get_id())
                                    catalyst_replacements[termx] = reactions
                                else:
                                    catalyst_replacements[termx] = [rule.get_id()]
                            if termx not in terms:
                                terms.append(termx)
                                mod_triplets.append([rxn, component, res])
                                if component in oppress_empty_residues.keys():
                                    ks = oppress_empty_residues[component]
                                    try: ks.remove(res)
                                    except: pass
                                    if ks != []:
                                        oppress_empty_residues[component] = ks
                                    else: oppress_empty_residues.pop(component,
                                                                     None)
                            make_lumped_contingencies(rule, modifier, rxn,
                                                      component, res)
                        elif res_l[1] is not None and res_r[1] is None:
                            try: rxn = mns[res_l[1].lower()]
                            except: rxn = 'Mod-'
                            rule.mod_type = rxn
                            res = res_l[0]
                            if res is None: res = residue
                            if '@' in component:
                                component = component.split('@')[0]
                            termx = '%s_%s_%s_[(%s)],%s,,,%s,%s,,%s,,,' % (modifier,
                                                                           rxn,
                                                                           component,
                                                                           res,
                                                                           modifier,
                                                                           rxn,
                                                                           component,
                                                                           res)
                            if '|' in modifier:
                                if termx in catalyst_replacements.keys():
                                    reactions = catalyst_replacements[termx]
                                    reactions.append(rule.get_id())
                                    catalyst_replacements[termx] = reactions
                                else:
                                    catalyst_replacements[termx] = [rule.get_id()]
                            if termx not in terms:
                                terms.append(termx)
                                mod_triplets.append([rxn, component, res])
                                if component in oppress_empty_residues.keys():
                                    ks = oppress_empty_residues[component]
                                    try: ks.remove(res)
                                    except: pass
                                    if ks != []:
                                        oppress_empty_residues[component] = ks
                                    else: oppress_empty_residues.pop(component,
                                                                     None)
                            make_lumped_contingencies(rule, modifier, rxn,
                                                      component, res)
                        else:
                            pass

        return terms

    def make_lumped_contingencies(rule, modifier, rxn, component, res):
        '''
        creates lumped contingencies for each modification;
        e.g. ABC -> ABC* would require !-contingencies AB and BC to execute *
        '''
        global contingencies
        global merge_instances
        global initialisations
        #global rxn2typeComp

        rxn2state = {'P+': 'P', 'Ub+': 'Ub', 'Palm+': 'Palm', 'Ace+': 'Ace',
                     'Pre+': 'pre', 'P-': '0', 'Ub-': '0', 'Palm-': '0',
                     'Ace-': '0', 'Pre-': '0'}
        instance_identifier = '<%s%s%s_%s>' % (modifier, rxn, component,
                                               rule.get_id())
        statement = '%s_%s_%s_[(%s)]' % (modifier, rxn, component, res)
        #rxn2typeComp[rule.get_id()] = [rxn,component]

        if statement in merge_instances.keys():
            if not type(merge_instances[statement]) == list:
                s = [merge_instances[statement]]
            else: s = merge_instances[statement]
            if len(s) > 1:
                instances = []
                for elem in merge_instances[statement]:
                    instances.append(elem)
            else: instances = [merge_instances[statement]]
            instances.append(instance_identifier)
            merge_instances[statement] = instances
        else: merge_instances[statement] = instance_identifier

        for component2 in rule.lhs:
            # THIS IS COMMENTED OUT FOR NOW SINCE IT SEEMED TO BE REDUNDANT
            # (SEE BUG#6 IN BITBUCKET)

            # required as the contingencies of this reaction are either
            # the ppis  and (NOT!) the statement itself
            # if len(set(component2.get_ppis())) > 1:
            #    for ppi in component2.get_ppis():
            #        contingencies.append('%s,1AND,%s' % (instance_identifier,
            #                                            ppi))
            # or a single ppi and the statement itself
            # el
            if len(set(component2.get_ppis())) == 1:
                contingencies.append('%s,AND,%s' % (instance_identifier,
                                                    component2.get_ppis()[0]))
            try:
                inits = initialisations[rule.get_id()]
                for i in inits:
                    if component2.get_class().lower() == 'complex':
                        equivalencies = create_equivalencies(component2)
                        contingencies.append('%s,AND,%s%s' % (instance_identifier,
                                                              i, equivalencies))
                    else:
                        contingencies.append('%s,AND,%s' % (instance_identifier,
                                                            i))
            except:
                pass

    def get_mods(element, mods=False):
        '''
        receives an entity and returns all modification residues of the
        topology
        '''
        mods = {}

        if element.get_class().lower() == 'complex':
            if element.get_definition():
                mods[element.get_name()] = element.get_definition().get_modification_residues()
            for entity in element.get_molecules():
                if entity.get_class().lower() == 'complex':
                    for em in entity.get_molecules():
                        more_mods = get_mods(em)
                        for mm in more_mods: mods[mm] = more_mods[mm]
                else:
                    try:
                        if entity.get_definition():
                            mods[entity.get_name()] = entity.get_definition().get_modification_residues()
                    except: pass
        elif element.get_kmer() == '2':
            mods[element.get_name() + '@' + str(element.get_index())] = element.get_definition().get_modification_residues()
            mods[element.get_name() + '@' + str(element.get_index() + 1)] = element.get_definition().get_modification_residues()
        else:
            if element.get_definition(): mods[element.get_name()] = element.get_definition().get_modification_residues()

        return mods

    def find_pairs(lhs, rhs):
        '''
        takes the left hand side and the right hand side and finds out which
        components are topologically identical, so we can then compare the
        modification residues
        '''
        lhs_topos = []
        rhs_topos = []
        pairs = []
        for elem in lhs: lhs_topos.append(make_topo(elem))
        for elem in rhs: rhs_topos.append(make_topo(elem))

        for i, ltopo in enumerate(lhs_topos):
            if ltopo in rhs_topos:
                for j, rtopo in enumerate(rhs_topos):
                    if ltopo == rtopo:
                        pairs.append([lhs[i], rhs[j]])

        return pairs

    def make_topo(element, topo=False):
        '''
        creates a characterising topology on ground of complex elements and
        names
        '''
        topo = ['<' + element.get_name() + '>']

        if element.get_class().lower() == 'complex':
            for entity in element.get_molecules():
                if entity.get_class().lower() == 'complex':
                    ss = make_topo(entity, topo)
                    for s in ss: topo.append(s)
                else:
                    try: topo.append(entity.get_protein_ref())
                    except: topo.append(entity.get_name())

        return sorted(topo)

    def create_lhsp(rule):
        '''
        creates the lhsP if required to unlump reactions
        '''
        if rule.rxn_type is False: return None
        elif rule.mod_type is False: return None
        else:
            if rule.synth != []:
                for elem in rule.lhs: rule.lhsp.append(elem)
                for elem in rule.synth:
                    blank_elem = copy.deepcopy(elem)
                    blank_elem.get_definition().clear_modification_residues()
                    rule.lhsp.append(blank_elem)
            elif rule.deg != []:
                for elem in rule.rhs: rule.lhsp.append(elem)
                for elem in rule.deg:
                    blank_elem = copy.deepcopy(elem)
                    blank_elem.get_definition().clear_modification_residues()
                    rule.lhsp.append(blank_elem)

    def determine_state_changes(rule):
        '''
        determines all the state changes between lhs, lhsP, and rhs
        '''
        global log_file
        global writelines
        global contingencies
        global catalyst_replacements
        global oppress_empty_residues

        bool_name = '%s_%s' % (rule.get_id(), '1')
        pls = {'phosphorylated': 'P+', 'acetylated': 'Ace+',
               'ubiquitinated': 'Ub+', 'palmytoylated': 'palm+',
               'prenylated': 'pre+', 'autophosphorylated': 'AP+'}
        mns = {'phosphorylated': 'P-', 'acetylated': 'Ace-',
               'ubiquitinated': 'Ub-', 'palmytoylated': 'palm-',
               'prenylated': 'pre-', 'autophosphorylated': 'AP-'}
        pairs1 = find_pairs(rule.lhs, rule.lhsp)
        pairs2 = find_pairs(rule.lhsp, rule.rhs)
        # all states that do not change between lhs and lhsp: context
        # of synth/deg reactions
        states_ds = []
        # all states that do not change between lhsp and rhs: context
        # of prod/con reactions
        states_pc = []

        substrate_names = [s.get_name() for s in rule.get_single_substrates()]
        product_names = [p.get_name() for p in rule.get_single_products()]
        degraded_names = []
        synth_names = []

        for s in substrate_names:
            if s not in product_names:
                degraded_names.append(s)
        for p in product_names:
            if p not in substrate_names:
                synth_names.append(p)

        for pair in pairs1:
            mods_l = get_mods(pair[0])        # modifications on lhs site
            mods_lp = get_mods(pair[1])        # modifications on lhs' site
            # save pair (incl. modifications; lets see what we will need here)
            if mods_l == mods_lp:
                states_ds.append(pair)

        for pair in pairs2:
            mods_lp = get_mods(pair[0])        # modifications on lhs' site
            mods_r = get_mods(pair[1])        # modifications on rhs site
            if mods_lp != mods_r:
                for k in mods_lp.keys():
                    if mods_lp[k] == {} or mods_r[k] == {}: continue
                    lp_k = mods_lp[k]
                    r_k = mods_r[k]
                    state_lp = lp_k.popitem()
                    state_r = r_k.popitem()
                    if state_lp[1][1] in mns.keys():
                        mode = mns[state_lp[1][1]]
                        if state_lp[0] != 'None': res = state_lp[0]
                        else: res = lp_k.popitem(0)
                    if state_r[1][1] in pls.keys():
                        mode = pls[state_r[1][1]]
                        if state_r[0] != 'None': res = state_r[0]
                        else: res = r_k.popitem(0)
                    else: continue
                    try: modifier = rule.get_modifiers()[0].get_name()
                    except: continue
                    if '@' in modifier: modifier = modifier.split('@')[0]
                    if '@' in k: k = k.split('@')[0]
                    row = '%s_%s_%s_[(%s)],%s,,,%s,%s,,%s,,,\n' % (modifier,
                                                                   mode, k, res,
                                                                   modifier,
                                                                   mode, k, res)
                    log_file += row
                    term = row
                    writelines.append(term)
                    if k in oppress_empty_residues.keys():
                        ks = oppress_empty_residues[k]
                        try: ks.remove(res)
                        except: pass
                        if ks != []: oppress_empty_residues[k] = ks
                        else: oppress_empty_residues.pop(k, None)
                    if '|' in modifier:
                        if term in catalyst_replacements.keys():
                            reactions = catalyst_replacements[term]
                            reactions.append(rule.get_id())
                            catalyst_replacements[term] = reactions
                        else: catalyst_replacements[term] = [rule.get_id()]
                    make_lumped_contingencies(rule, modifier, mode, k, res)
            else:
                # save pair (incl. modifications; lets see what we need here)
                states_pc.append(pair)

        # generate BOOLs for synthesis/degradation
        # if synth_names != [] and diff_states != []:
        #     try: modifier = rule.get_modifiers()[0].get_name()
        #     except: modifier = 'UnknownComponent'
        #     for sn in synth_names:
        #         writelines.append('%s_[DomainA]_Syn_%s_[DomainB],%s,[Domain'\
        #                            'A],[ResidueA],Syn,%s,[DomainB],[DomainB'\
        #                            '],,,' % (modifier, sn, modifier, sn))
        #         contingencies.append('%s,AND,%s_[DomainA]_Syn_%s_'\
        #                              '[DomainB]' % (modifier, sn))
        #         print('%s,AND,%s_[DomainA]_Syn_%s_[DomainB]' % (bool_name,
        #                                                          modifier, sn))
        #     for dn in degraded_names:
        #         writelines.append('%s_[DomainA]_Deg_%s_[DomainB],%s,[Domain'\
        #                           'A],[ResidueA],Deg,%s,[DomainB],[DomainB]'\
        #                           ',,,' % (modifier, dn, modifier, dn))
        #         contingencies.append('%s,AND,%s_[DomainA]_Deg_%s_[DomainB'\
        #                              ']' % (modifier, dn))
        #         print('%s,AND,%s_[DomainA]_Deg_%s_[DomainB]' % (bool_name,
        #                                                         modifier,
        #                                                         dn))

    def gather_res(mi):
        '''
        collect residues for multiple modification reactions
        '''
        residues = []
        for elem in mi:
            residues.append(re.search('(.*)\\[\\((.*)\\)\\]', elem).group(2))
            try: residues.append(re.search('(.*)\\[\\((.*)\\)\\]',
                                           elem).group(2))
            except: pass
        return set(residues)

    def analyse_merge_instances():
        '''
        analyse whether certain reaction instances are duplicate and can be
        merged
        '''
        global merge_instances
        global contingencies

        for elem in merge_instances.keys():
            makebool = copy.deepcopy(elem)
            for sc in ['_', '[', ']']: makebool = makebool.replace(sc, '')
            contingencies.append('%s,!,<%s>' % (elem, makebool))
            mi = merge_instances[elem]
            if not type(mi) == list: mi = [mi]
            if len(set(mi)) > 1:
                for m in mi:
                    contingencies.append('<%s>,OR,%s' % (makebool, m))
            else:
                contingencies.append('<%s>,AND,%s' % (makebool, mi[0]))

    def create_equivalencies(component):
        '''
        creates the equivalencies for the complex structure contingencies,
        e.g. #Ste11@1=Ste11@4
        '''
        global log_file
        equivalencies = []

        for i, element in enumerate(component.get_molecules()):
            if element.get_class().lower() == 'complex':
                equivalencies_iter = create_equivalencies(element)
                if equivalencies_iter not in equivalencies:
                    equivalencies.append(equivalencies_iter)
                #equivalencies.append(create_equivalencies(element))
            else:
                if not element.get_index_int(): element.set_index_int(i+1)
                try:
                    equivalency = '#%s@%s=%s@%s' % (element.get_name(),
                                                           element.get_index(),
                                                           element.get_name(),
                                                           element.get_index_int())
                    if equivalency not in equivalencies:
                        equivalencies.append(equivalency)


                except:
                    log_file += 'Warning: Couldnt create an equivalency for %s in complex'\
                                '%s.\n' % (element.get_name(), component.get_name())

        return ''.join(equivalencies)

    def find_inherited_equivalencies(rule):
        '''
        create all inherited equivalencies of all components of a rule
        '''
        global log_file
        equivalencies = []

        for i, element in enumerate(rule.get_left_hand_side()):
            if element.get_class().lower() == 'complex':
                for i, element2 in enumerate(element.get_molecules()):
                    if element2.get_class().lower() == 'complex':
                        equivalencies_iter = create_equivalencies(element2)
                        if equivalencies_iter not in equivalencies:
                            equivalencies.append(equivalencies_iter)
                    else:
                        if not element2.get_index_int(): element2.set_index_int(i+1)
                        try:
                            equivalency = '#%s@%s=%s@%s' % (element2.get_name(),
                                                                   element2.get_index(),
                                                                   element2.get_name(),
                                                                   element2.get_index_int())
                            if equivalency not in equivalencies:
                                equivalencies.append(equivalency)
                        except:
                            log_file += 'Warning: Couldnt create an equivalency for %s in complex'\
                                        '%s.\n' % (element2.get_name(), element2.get_name())
            else:
                if not element.get_index_int(): element.set_index_int(i+1)
                try:
                    equivalency = '#%s@%s=%s@%s' % (element.get_name(),
                                                    element.get_index(),
                                                    element.get_name(),
                                                    element.get_index_int())
                    if equivalency not in equivalencies:
                        equivalencies.append(equivalency)
                except:
                    log_file += 'Warning: Couldnt create an equivalency for %s in complex'\
                                '%s.\n' % (element.get_name(),
                                           element.get_name())

        return ''.join(equivalencies)

    
    def initialise_states(rule):
        '''
        if a reaction has a component with a state (modified or unmodified),
        this needs to be initialised as a contingency
        '''
        global contingencies
        global initialisations
        global mod_triplets_new
        global oppress_empty_residues

        long2short = {'phosphorylated': 'P', 'acetylated': 'Ace',
                      'ubiquitinated': 'Ub', 'palmytoylated': 'palm',
                      'prenylated': 'pre', None: '0'}

        for component in rule.get_left_hand_side():
            component_mods = get_mods(component)
            if not component.get_class().lower() != 'complex':
                equivalencies = create_equivalencies(component)
                contingencies.append('<%s_%s>,AND,<%s>%s' % (component.get_name(),
                                                             rule.get_id(),
                                                             component.get_name(),
                                                             equivalencies))
                if rule.get_id() in initialisations.keys():
                    inits = initialisations[rule.get_id()]
                    inits.append('<%s_%s>' % (component.get_name(),
                                              rule.get_id()))
                    initialisations[rule.get_id()] = inits
                else:
                    initialisations[rule.get_id()] = ['<%s_%s>' % (component.get_name(),
                                                                   rule.get_id())]

            if component.get_definition() and component.get_definition().get_modification_residues() == {}:
                continue
            for k in component_mods.keys():
                k_obj = None
                if component_mods[k] == {}: continue
                single_dict = component_mods[k]
                for obj in find_components([component]):
                    if obj.get_name() == k:
                        k_obj = obj
                        continue
                    if '@' in obj.get_name() and obj.get_name().split('@')[0] == k:
                        k_obj = obj
                        continue
                    if '@' in k and obj.get_name() == k.split('@')[0]:
                        k_obj = obj

                if not k_obj:
                    log_file += 'Warning: There is a problem with element'\
                                '%s in rule %s.\n' % (k, rule.get_id())
                    continue

                for l in single_dict.keys():
                    if single_dict[l][0]: residue = single_dict[l][0]
                    else: residue = l
                    residue_content = single_dict[l]
                    if k_obj.get_class().lower() == 'complex':
                        contingency = '<%s_%s>,AND,%s_[(%s)]-{%s}' % (component.get_name(),
                                                                      rule.get_id(), k,
                                                                      residue,
                                                                      long2short[residue_content[1]])
                    else:
                        if '@' not in k:
                            contingency = '<%s_%s>,AND,%s@%s_[(%s)]-{%s}' % (component.get_name(),
                                                                             rule.get_id(), k,
                                                                             k_obj.get_index(),
                                                                             residue,
                                                                             long2short[residue_content[1]])
                        else:
                            contingency = '<%s_%s>,AND,%s_[(%s)]-{%s}' % (component.get_name(),
                                                                          rule.get_id(), k,
                                                                          residue,
                                                                          long2short[residue_content[1]])
                        if long2short[residue_content[1]] == '0':
                            if '@' in k: p = k.split('@')[0]
                            else: p = k
                            if p in oppress_empty_residues.keys():
                                ks = oppress_empty_residues[p]
                                if residue not in ks: ks.append(residue)
                                oppress_empty_residues[p] = ks
                            else: oppress_empty_residues[p] = [residue]
                    contingencies.append(contingency)
                    if not long2short[residue_content[1]] == '0':
                        try: k = k.split('@')[0]
                        except: pass
                        mod_triplets_new.append([long2short[residue_content[1]] + '+',
                                                 k, residue])
                    if rule.get_id() in initialisations.keys():
                        inits = initialisations[rule.get_id()]
                        inits.append('<%s_%s>' % (component.get_name(),
                                                  rule.get_id()))
                        initialisations[rule.get_id()] = inits
                    else:
                        initialisations[rule.get_id()] = ['<%s_%s>' % (component.get_name(),
                                                                       rule.get_id())]
    def complement_reactions(ppi_duplets, ppi_duplets_new,
                             mod_triplets, mod_triplets_new):
        '''
        complement the reaction list according to reactions that are not
        explicitly listed in the model, but which are required for the rxncon
        processing; e.g. a ppi of two components which are always connected in
        the model, but the interaction is no official reaction
        '''
        global writelines
        
        # start with the ppis:
        for ppi_duplet in ppi_duplets:
            line = '%s_[%s]_ppi+_%s_[%s],%s,%s,,ppi+,%s,%s,,,,' % (ppi_duplet[0],
                                                                   ppi_duplet[1],
                                                                   ppi_duplet[1],
                                                                   ppi_duplet[0],
                                                                   ppi_duplet[0],
                                                                   ppi_duplet[1],
                                                                   ppi_duplet[1],
                                                                   ppi_duplet[0])
            for k in ['GTP', 'GDP', 'GMP', 'ATP', 'ADP', 'AMP']:
                if k in ppi_duplet: line = line.replace('ppi', 'i')
            writelines.append(line)

        # then the modifications
        for mod_triplet in mod_triplets_new:
            if mod_triplet not in mod_triplets:
                line = '%s_%s_%s_[%s],%s,,,%s,%s,%s,,,,' % ('UnknownModifier',
                                                            mod_triplet[0],
                                                            mod_triplet[1],
                                                            mod_triplet[2],
                                                            'UnknownModifier',
                                                            mod_triplet[0],
                                                            mod_triplet[1],
                                                            mod_triplet[2])
            writelines.append(line)

    def construct_ppi_contingencies(rule, ppi1, ppi2, ppi_type):
        '''
        if the reaction is an interaction, we need to define the corresponding
        contingencies for it.
        '''
        global contingencies

        contingencies.append('%s_%s_%s,!,<%s%s%s>' % (ppi1, ppi_type, ppi2,
                                                    ppi1, ppi_type, ppi2))


        # find equivalencies for everything in rule.getlefthand
        
        equivalencies = find_inherited_equivalencies(rule)
        contingencies.append('<%s%s%s>,AND,<%s%s%s_%s>%s' % (ppi1, ppi_type, ppi2,
                                                              ppi1, ppi_type, ppi2,
                                                              rule.get_id(),
                                                              equivalencies))

        for component in rule.get_left_hand_side():
            if component.get_class().lower() == 'complex':
                equivalencies = create_equivalencies(component)
            else: equivalencies = find_inherited_equivalencies(rule)

            # find equivalencies for everything under/above rule.getlefthand
            
            contingencies.append('<%s%s%s_%s>,AND,<%s_%s>%s' % (ppi1, ppi_type, ppi2,
                                                                 rule.get_id(),
                                                                 component.get_name(),
                                                                 rule.get_id(),
                                                                 equivalencies))

            if component.get_class().lower() != 'complex':
                contingencies.append('<%s_%s>,AND,%s@%s' % (component.get_name(),
                                                            rule.get_id(),
                                                            component.get_name(),
                                                            component.get_index()))
                

    ########################
    # End of subreactions ##
    ########################
    global log_file
    global writelines
    global events
    global contingencies
    global merge_instances
    global r2terms
    global catalyst_replacements
    global reaction2singlecat
    global singlecat2residue
    global initialisations
    global ppi_duplets
    global mod_triplets
    global mod_triplets_new
    global oppress_empty_residues
    log_file = ''
    #global rxn2typeComp
    ppi_duplets_new = []
    mod_triplets = []
    mod_triplets_new = []
    writelines = []
    events = {}
    merge_instances = {}
    r2terms = {}
    initialisations = {}
    catalyst_replacements = {}
    reaction2singlecat = {}
    singlecat2residue = {}
    oppress_empty_residues = {}
    #rxn2typeComp = {}

    # 1: prepare the conversion by adding knowledge about complex topology
    model_object, ppi_duplets, contingencies = prepare(file_name)

    # 2: find out reaction types by iterating over the reactions
    #    and comparing (i) the components, (ii) their ppis, and (iii) their
    #    modification states
    for rule in model_object.get_rules():
        # 2.1 see if there are differences in the components
        lhs = rule.get_left_hand_side()
        rhs = rule.get_right_hand_side()

        try: modifier = rule.get_modifiers()[0].get_name()
        except: modifier = 'No Modifier'
        log_file += '\n############\n%s\nReaction:' % (rule.get_id()) + \
                    str([su.get_name() for su in lhs]) + '->' + \
                    str([sa.get_name() for sa in rhs]) + '#' + modifier + '\n'

        for elem in lhs:
            c = find_components([elem])
            for d in c: rule.add_single_substrate(d)
        for elem in rhs:
            c = find_components([elem])
            for d in c: rule.add_single_product(d)
        rule.sort_entities()
        log_file += 'Single components:' + \
                    str(sorted([k.get_name() + '@' + str(k.get_index()) for k in rule.get_single_substrates()])) + \
                    '->' + \
                    str(sorted([p.get_name() + '@' + str(p.get_index()) for p in rule.get_single_products()])) + \
                    '\n'
        # (2.1a check for phenotypes)
        # eve = pheno_check(rule, lhs, rhs)

        # (2.1b check for transc, transl, transp, etc.)
        trans = trans_check(rule)
        if trans is not None:
            if type(trans) == list:
                for elem in trans: writelines.append(elem)
                continue
            else:
                writelines.append(trans)
                continue

        initialise_states(rule)

        delta_c(rule)
        log_file += 'Rxn type: %s\n' % (rule.rxn_type)

        # 2.2 see if there are differences in the components' ppis
        lhs_ppis = []
        rhs_ppis = []
        for comp in lhs:
            if comp.get_class().lower() != 'complex':
                if comp.get_kmer() != '2':
                    continue
            pps = comp.get_ppis()
            for c in pps: lhs_ppis.append(c)
        for comp in rhs:
            if comp.get_class().lower() != 'complex':
                if comp.get_kmer() != '2':
                    continue
            pps = comp.get_ppis()
            for c in pps: rhs_ppis.append(c)
            
        [lost_ppis, new_ppis] = delta_ppi(lhs_ppis, rhs_ppis)
        for lost_ppi in lost_ppis:
            try:
                ppi1 = re.search('(.*)--', lost_ppi).group(1)
                ppi2 = re.search('--(.*)', lost_ppi).group(1)
                if '@' in ppi1: ppi1 = ppi1.split('@')[0]
                if '@' in ppi2: ppi2 = ppi2.split('@')[0]
                if ppi1 in model_object.get_simple_molecules():
                    row = '%s_[%s]_i-_%s_[%s],%s,%s,,i-,%s,%s,,,,' % (ppi2,
                                                                      ppi1,
                                                                      ppi1,
                                                                      ppi2,
                                                                      ppi2,
                                                                      ppi1,
                                                                      ppi1,
                                                                      ppi2)
                    log_file += row + '\n'
                    writelines.append(row)
                    construct_ppi_contingencies(rule, ppi1, ppi2, ppi_type='i-')
                elif ppi2 in model_object.get_simple_molecules():
                    row = '%s_[%s]_i-_%s_[%s],%s,%s,,i-,%s,%s,,,,' % (ppi1,
                                                                      ppi2,
                                                                      ppi2,
                                                                      ppi1,
                                                                      ppi1,
                                                                      ppi2,
                                                                      ppi2,
                                                                      ppi1)
                    log_file += row + '\n'
                    writelines.append(row)
                    construct_ppi_contingencies(rule, ppi2, ppi1, ppi_type='i-')
                else:
                    row = '%s_[%s]_ppi-_%s_[%s],%s,%s,,ppi-,%s,%s,,,,' % (ppi1,
                                                                          ppi2,
                                                                          ppi2,
                                                                          ppi1,
                                                                          ppi1,
                                                                          ppi2,
                                                                          ppi2,
                                                                          ppi1)
                    log_file += row + '\n'
                    writelines.append(row)
                    construct_ppi_contingencies(rule, ppi1, ppi2, ppi_type='ppi-')
            except: log_file += 'Warning: Could not articulate lost ppi %s' % (lost_ppi)

        for new_ppi in new_ppis:
            try:
                ppi1 = re.search('(.*)--', new_ppi).group(1)
                ppi2 = re.search('--(.*)', new_ppi).group(1)
                if '@' in ppi1: ppi1 = ppi1.split('@')[0]
                if '@' in ppi2: ppi2 = ppi2.split('@')[0]
                if ppi1 in model_object.get_simple_molecules():
                    row = '%s_[%s]_i+_%s_[%s],%s,%s,,i+,%s,%s,,,,' % (ppi2,
                                                                      ppi1,
                                                                      ppi1,
                                                                      ppi2,
                                                                      ppi2,
                                                                      ppi1,
                                                                      ppi1,
                                                                      ppi2)
                    log_file += row + '\n'
                    writelines.append(row)
                    ppi_duplets_new.append([ppi2, ppi1])
                    construct_ppi_contingencies(rule, ppi2, ppi1, ppi_type='i+')
                elif ppi2 in model_object.get_simple_molecules():
                    row = '%s_[%s]_i+_%s_[%s],%s,%s,,i+,%s,%s,,,,' % (ppi1,
                                                                      ppi2,
                                                                      ppi2,
                                                                      ppi1,
                                                                      ppi1,
                                                                      ppi2,
                                                                      ppi2,
                                                                      ppi1)
                    log_file += row + '\n'
                    writelines.append(row)
                    ppi_duplets_new.append([ppi1, ppi2])
                    construct_ppi_contingencies(rule, ppi1, ppi2, ppi_type='i+')
                else:
                    row = '%s_[%s]_ppi+_%s_[%s],%s,%s,,ppi+,%s,%s,,,,' % (ppi1,
                                                                          ppi2,
                                                                          ppi2,
                                                                          ppi1,
                                                                          ppi1,
                                                                          ppi2,
                                                                          ppi2,
                                                                          ppi1)
                    log_file += row + '\n'
                    writelines.append(row)
                    ppi_duplets_new.append([ppi1, ppi2])
                    construct_ppi_contingencies(rule, ppi1, ppi2, ppi_type='ppi+')
            except: log_file += 'Warning: Could not articulate new ppi %s' % (new_ppi)


        # 2.3 see if there are differences in components' modification states
        lhs_mods = []
        rhs_mods = []
        for comp in lhs:
            try:
                comp_def = comp.get_definition()
                mods = comp_def.get_modification_residues()
                for c in pps: lhs_ppis.append(c)
            except: pass
        for comp in rhs:
            try:
                comp_def = comp.get_definition()
                mods = comp_def.get_modification_residues()
                for c in pps: lhs_ppis.append(c)
            except: pass

        terms = delta_m(rule, lhs, rhs)
        if terms != []:
            for termx in terms:
                writelines.append(termx)
        log_file += 'Modification type: %s\n' % (rule.mod_type)

        create_lhsp(rule)
        log_file +='LHSp:' + str([k.get_name() for k in rule.lhsp]) + '\n'
        if rule.lhsp != []: determine_state_changes(rule)

    analyse_merge_instances()
    # add the additional writelines for reactions that do not appear in the
    # model but are required
    complement_reactions(ppi_duplets, ppi_duplets_new, mod_triplets,
                         mod_triplets_new)

    new_wl = get_catalysts_by_user(writelines)

    for rule in model_object.get_rules():
        if rule.get_id() in reaction2singlecat.keys():
            rule.set_modifier_name(reaction2singlecat[rule.get_id()])

    # if events != {}:
    new_wl.append('\n\n!!SBtab TableType="contingencies" TableName="Contingency'\
                  'List"\n!ID,!Target,!Contingency,!Modifier')

    for e in events:
        for s in events[e]:
            try:
                if len(events[e]) > 2: contingencies.append(''.join(s))
                else: contingencies.append(''.join(s.replace('OR', 'AND')))
            except: log_file += 'Error writing event %s\n' % (e)

    outwriter(new_wl, contingencies)

    return model_object


def get_catalysts_by_user(writelines):
    '''
    output and input required by user to determine the correct catalyst of a
    reaction what we have is something like this: "[A|B|C]_P+_D", what we need
    is "A_P+_D"
    '''
    global reaction2singlecat
    global singlecat2residue
    global catalyst_replacements

    ff = re.search('(.*)/(.*)', file_name[:-4].lower()).group(2)
    filename = 'files/%s_catalysts.csv' % (ff)
    f = open(filename, 'w')
    f.write('!!SBtab TableType="rxnconReactionList" TableName="ReactionList"\n'\
            '!ID,!UID:Reaction,!ComponentA:Name,!ComponentA:Domain,!ComponentA'\
            ':Residue,!Reaction,!ComponentB:Name,!ComponentB:Domain,!Component'\
            'B:Residue,!Quality,!Literature:Identifiers:pubmed,!Comment\n')
    dismiss = []
    old = []
    counter = 1
    for i, line in enumerate(writelines):
        if line not in dismiss:
            f.write(str(counter) + ',' + line + '\n')
            old.append(line)
            dismiss.append(line)
            counter += 1
    f.close()

    print('The file catalysts.csv has been written to the files directory. P'\
          'lease open the file, choose a catalyst for each ambiguous reaction,'\
          'and save it. Please do not change the order of reactions or add/rem'\
          'ove rows. If you"re done, please hit enter and this script will pro'\
          'ceed.')
    input('')

    get_file = False
    while get_file is False:
        try:
            catalysts_file = open(filename[:-4] + '_re.csv', 'r').read()
            get_file = True
        except:
            print('The catalysts file could not be read. Please make sure the'\
                  'file name and location have not changed and then hit enter'\
                  'again.')
            input('')

    new_wl = []
    auto_p = False
    for row in catalysts_file.split('\n'):
        if '|' in row:
            print('There is an error in row %s. Please choose the correct cat'\
                  'alyst in the response file and run the program'\
                  'again.' % (row))
            input('')
        else: new_wl.append(row)
        if 'AP+' in row:
            row = row.replace('AP+', 'P+')
            auto_p = True

        # see output:
        # the rest == rest2 comparison is too ambiguous; it fits for all
        # reX_Y scenarios and thus also overwrites all of them in the end.
        splitrow = ','.join(row.split(',')[1:])
        if splitrow in old or splitrow.startswith('!') or splitrow == '' or splitrow.startswith(','): continue
        else:
            splitrow = splitrow.split(',')
            try:
                residue = re.search(r'(\((.*?)\))', splitrow[0]).group(0).strip('()')
            except: residue = None
            catalyst = splitrow[1]
            rest = splitrow[2:8]
            # here in the rest: if rest[1] == rest[5]
            # --> this is an autophosphorylation --> rest[4] = 'AP+'
            for cr in catalyst_replacements.keys():
                catalyst_string = cr.split(',')[1]
                rest2 = cr.split(',')[2:8]
                if rest == rest2 and catalyst in catalyst_string:
                    for elem in catalyst_replacements[cr]:
                        reaction2singlecat[elem] = catalyst
                        if catalyst in singlecat2residue.keys():
                            residues = singlecat2residue[catalyst]
                            if residue not in residues:
                                residues.append(residue)
                            singlecat2residue[catalyst] = residues
                        else: singlecat2residue[catalyst] = [residue]

    return new_wl


def check_boolean_inits(outputs):
    '''
    check whether all boolean expressions (<..>) in the contingencies are
    also defined as targets.
    currently, we remove expressions, which are not initialised alongside a
    warning message; the check function iterates over the list until there
    are no more uninitialised expressions.
    '''
    global log_file
    error = True

    while error:
        targets = []
        modifiers = []
        mod2line = {}
        tar2line = {}
        error = False
        remove_lines = []
        for op in outputs:
            splits = op.split(',')
            if splits[0].startswith('<'):
                target = splits[0]
                targets.append(target)
                if target in tar2line:
                    t = tar2line[target]
                    t.append(op)
                    tar2line[target] = t
                else: tar2line[target] = [op]
                
            if splits[2].startswith('<'):
                if '#' in splits[2]: s = splits[2].split('#')[0].rstrip('\n')
                else: s = splits[2].rstrip('\n')
                modifiers.append(s)
                if s in mod2line:
                    ss = mod2line[s]
                    ss.append(op)
                    mod2line[s] = ss
                else:
                    mod2line[s] = [op]

        for m in set(modifiers):
            if m not in set(targets):
                log_file += 'Warning: The modifier %s has not been initialised. Please '\
                            'initialise or remove line %s\n' % (m, mod2line[m])
                remove_lines.append(mod2line[m])
                mod2line.pop(m)
                error = True

        for line in remove_lines:
            for entry in line:
                outputs.remove(entry)
                log_file += 'Line %s has been removed.\n' % (entry)
    
    return (outputs, targets, modifiers, mod2line, tar2line)


def flatten_contingencies(outputs, targets, modifiers, mod2line, tar2line):
    '''
    contingencies with only one child can be flattened; and while we're at
    it, we can try and inherit the equivalencies
    '''
    flattened = True
    remove_pairs = []
    
    while flattened:
        flattened = False
        for target in targets:
            if modifiers.count(target) == 1 and targets.count(target) == 1:
                splits = tar2line[target][0].split(',')
                new_target = splits[2]
              
                if splits[2].startswith('<'):
                    for modline in mod2line[target]:
                        if '#' in modline: modline = modline.split('#')[0]
                        new_mod_line = modline.replace(target,new_target)
                        if new_mod_line not in outputs:
                            outputs.append(new_mod_line)
                        try:
                            outputs.remove(modline)
                        except: pass
                        
                        if new_target in mod2line:
                            lines = mod2line[new_target]
                            lines.append(new_mod_line)
                            mod2line[new_target] = lines
                        else: mod2line[new_target] = [new_mod_line]
                    mod2line.pop(target)
                    
                    for target_line in tar2line[target]:
                        if new_target in tar2line:
                            lines = tar2line[new_target]
                            lines.append(new_mod_line)
                            tar2line[new_target] = lines
                        else: tar2line[new_target] = [new_mod_line]
                    tar2line.pop(target)
                    modifiers.remove(target)
                    flattened = True

    return outputs
    


def correct_contingency_lines(contingencies, dismiss):
    '''
    write the contingencies and dont forget to replace certain
    characters and ambiguities
    '''
    global reaction2singlecat
    global singlecat2residue
    #global rxn2typeComp
    
    rmv_special_chars = {'P\\+': 'Pplus', 'Ub\\+': 'Ubplus',
                         'Palm\\+': 'Palmplus', 'Ace\\+': 'Aceplus',
                         'Pre\\+': 'preplus', 'P\\-': 'Pminus',
                         'Ub\\-': 'Ubminus', 'Palm\\-': 'Palmminus',
                         'Ace\\-': 'Aceminus', 'Pre\\-': 'Preminus',
                         'ppi\\+': 'ppiplus', 'ppi\\-': 'ppiminus',
                         'i\\+': 'iplus', 'i\\-': 'iminus',
    }
    modifications = ['P+','Ub+','Palm+','Pre+']
    
    #counter = 1
    outputs = []
    for j, line in enumerate(contingencies):
        empty_residue = False
        if line not in dismiss:
            # 1: cutout empty residues
            for component in oppress_empty_residues.keys():
                for residue in oppress_empty_residues[component]:
                    if component in line and '_[(' + residue + ')]-{0}' in line:
                        empty_residue = True
            if empty_residue:  continue
            '''
            # 1b: cutout source states from modification reactions
            for typeCompKey in rxn2typeComp.keys():
                if typeCompKey+'>' in line and  ')]-{0}' in line:
                    tC = rxn2typeComp[typeCompKey]
                    #print(tC)
                    #input()
                    if tC[0] in modifications:
                        print(line)
                        print(tC)
                        input()
            '''
            # 2: replace ambiguous catalysts with correct catalysts for the
            #    different cases
            if '|' in line:
                for r_id in reaction2singlecat.keys():
                    if r_id in line:    # case: reaction ID appears in line
                        new_catalyst = reaction2singlecat[r_id]
                        old_catalyst = re.search(r'(\[(.*?)\])', line)
                        try:
                            line = line.replace(old_catalyst.group(1),
                                                new_catalyst)
                            line = line.replace(old_catalyst.group(1).strip('[]'),
                                                new_catalyst)
                            if '_TRSC_' in line:
                                line = line.replace('_' + r_id, '')
                            break
                        except: pass
                    # case: reaction ID with underscoreDigit appears in line
                    if r_id.split('_')[0] in line:
                        try:
                            r_id_iter = re.search(r_id.split('_')[0] + '_.',
                                                  line).group(0)
                        except: r_id_iter = r_id
                        try:
                            try: new_catalyst = reaction2singlecat[r_id]
                            except:
                                new_catalyst = reaction2singlecat[r_id.split('_')[0]]
                            old_catalyst = re.search(r'(\[(.*?)\])', line)
                            try:
                                if new_catalyst in old_catalyst:
                                    line = line.replace(old_catalyst.group(1),
                                                        new_catalyst)
                                    line = line.replace(old_catalyst.group(1).strip('[]'),
                                                        new_catalyst)
                                    if '_TRSC_' in line:
                                        line = line.replace('_' + r_id_iter, '')
                                    break
                            except: pass
                        except: pass

                # in transcriptions the reaction identifier is required, but
                # needs to be removed in the end
                if '_TRSC_' in line:
                    try:
                        r_id = re.search('_re.*_.', line).group(0)
                        line = line.replace(r_id, '')
                    except: pass

                for cat in reaction2singlecat.values():
                    if cat in line:
                        breaker = False
                        old_catalyst = re.search(r'(\[(.*?)\])', line)
                        try:
                            residues = singlecat2residue[cat]
                            for r in residues:
                                if r in line:
                                    line = line.replace(old_catalyst.group(1),
                                                        cat)
                                    line = line.replace(old_catalyst.group(1).strip('[]'),
                                                        cat)
                                    breaker = True
                                    break
                            if breaker: break
                        except: pass

            # 3: remove special characters if the appear in
            #    Boolean statements <...>
            for sc in rmv_special_chars.keys():
                while re.search('<.*(%s).*>' % (sc), line):
                    terms = re.search('(.*)(<.*)(%s)(.*>)(.*)' % (sc), line)
                    ss = sc.replace('\\', '')
                    line = terms.group(1) + terms.group(2) + terms.group(3).replace(ss, rmv_special_chars[sc]) + terms.group(4) + terms.group(5)
            #nl = str(counter) + ',' + ''.join(line) + '\n'
            nl = ''.join(line)
            if line not in dismiss:
                if '|' in nl: nl = nl.replace('|', '')
                outputs.append(nl)
                dismiss.append(line)
                #counter += 1

    return outputs

            
def outwriter(writelines, contingencies):
    '''
    writes to hard disk
    '''
    ff = re.search('(.*)/(.*)', file_name[:-4].lower()).group(2)
    filename = 'files/%s_model.csv' % (ff)
    f = open(filename, 'w')
    dismiss = []
    
    # first write the reactions
    for i, line in enumerate(writelines):
        if line not in dismiss:
            f.write(line + '\n')
            dismiss.append(line)

    outputs = correct_contingency_lines(contingencies, dismiss)
    
    (outputs_init, targets, modifiers,
     mod2line, tar2line) = check_boolean_inits(outputs)
   
    outputs_flattened = flatten_contingencies(outputs_init, targets,
                                              modifiers, mod2line, tar2line)
    
    (outputs_init, targets, modifiers,
     mod2line, tar2line) = check_boolean_inits(outputs_flattened)

    doubles = []
    for i,op in enumerate(sorted(outputs_init)):
        if not op in doubles:
            f.write(str(i+1) + ',' + op + '\n')
            doubles.append(op)
    f.close()

    ll = re.search('(.*)/(.*)', file_name[:-4].lower()).group(2)
    filename = 'files/%s_logfile.csv' % (ll)
    l = open(filename,'w')
    for line in log_file:
        l.write(line)
    l.close()
    

def prepare(file_name: str):
    '''
    imports the SBML model, makes it a Python object, and adds
    topology knowledge
    '''
    (model_object, ppi_duplets, contingencies) = prepareConversion.prepareModel(file_name)
    return model_object, ppi_duplets, contingencies


if __name__ == '__main__':

    try: sys.argv[1]
    except: sys.exit()

    file_name = sys.argv[1]
    model_object = convert(file_name)
    print('')
