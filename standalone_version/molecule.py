#!/usr/bin/env python


class MoleculeDef:

    def __init__(self, ident: str, name: str):
        self.ident = ident
        self.name = name
        self.modification_residues = {}
        self.bindings = {}

    def set_notes(self, notes: str):
        self.notes = notes

    def set_complex(self, complex_name: str):
        self.parent_complex = complex_name

    def set_species_identity(self, identity: str):
        self.species_identity = identity

    def set_protein_reference(self, reference: str):
        self.protein_reference = reference

    def add_modification_residue(self, residue: str, state: str,
                                 name: str=False, angle=False, side=False):
        self.modification_residues[residue] = (name, state, angle, side)

    def get_modification_residues(self):
        return self.modification_residues

    def clear_modification_residues(self):
        for mr in self.modification_residues.keys():
            attributes = self.modification_residues[mr]
            new_attributes = (attributes[0], None, attributes[2],
                              attributes[3])
            self.modification_residues[mr] = new_attributes

    def set_mod_state(self, residue: str, state: str):
        mod_residue = self.modification_residues[residue]
        new_mod_residue = (mod_residue[0], state, mod_residue[2],
                           mod_residue[3])
        self.modification_residues[residue] = new_mod_residue

    def get_modification_residue(self, query: str):
        try: return self.modification_residues[query]
        except: return None

    def set_binding_domain(self, partner: str, domain: str):
        self.bindings[partner] = domain

    def get_binding_domains(self):
        return self.bindings


class ComplexDef:

    def __init__(self, ident: str, name: str):
        self.ident = ident
        self.name = name
        self.modification_residues = {}
        self.protein_reference = None
        self.bindings = {}
        self.species_identity = None

    def set_notes(self, notes: str):
        self.notes = notes

    def get_notes(self):
        return self.notes

    def set_complex(self, complex_name: str):
        self.parent_complex = complex_name

    def set_species_identity(self, identity: str):
        self.species_identity = identity

    def get_species_identity(self):
        return self.species_identity

    def set_protein_reference(self, reference: str):
        self.protein_reference = reference

    def get_protein_reference(self):
        return self.protein_reference

    def add_modification_residue(self, residue: str, state: str, name=False,
                                 angle=False, side=False):
        self.modification_residues[residue] = (name, state, angle, side)

    def get_modification_residues(self):
        return self.modification_residues

    def set_mod_state(self, residue: str, state: str):
        mod_residue = self.modification_residues[residue]
        new_mod_residue = (mod_residue[0], state, mod_residue[2],
                           mod_residue[3])
        self.modification_residues[residue] = new_mod_residue

    def get_modification_residue(self, query: str):
        try: return self.modification_residues[query]
        except: return None

    def set_binding_domain(self, partner: str, domain: str):
        self.bindings[partner] = domain

    def get_binding_domains(self):
        return self.bindings


class Molecule:

    def __init__(self, ident: str, name: str):
        self.ident = ident
        self.name = name
        self.cat_reactions = []
        self.stoichiometry = 1.
        self.molecule_def = False
        self.sclass = None
        self.kmer = 1
        self.ppis = []
        self.compartment = None
        self.position = False
        self.index = None
        self.index_int = None

    def set_index(self, index):
        self.index = index

    def get_index(self):
        return self.index

    def set_index_int(self, index):
        self.index_int = index

    def get_index_int(self):
        return self.index_int

    def get_id(self):
        return self.ident

    def get_name(self):
        return self.name

    def set_name(self, name):
        self.name = name

    def set_compartment(self, compartment):
        self.compartment = compartment

    def get_compartment(self):
        return self.compartment

    def set_position(self, position):
        self.position = position

    def get_position(self):
        return self.position

    def set_class(self, sclass: str):
        self.sclass = sclass

    def get_class(self):
        return self.sclass

    def set_name_cd(self, name_cd: str):
        self.name_cd = name_cd

    def get_name_cd(self):
        return self.name_cd

    def set_protein_ref(self, protein_ref: str):
        self.protein_ref = protein_ref

    def get_protein_ref(self):
        return self.protein_ref

    def set_molecule_def(self, MoleculeDef: MoleculeDef):
        self.molecule_def = MoleculeDef

    def get_definition(self):
        return self.molecule_def

    def set_cat_reactions(self, cat_reactions):
        self.cat_reactions.append(cat_reactions)

    def get_cat_reactions(self):
        return self.cat_reactions

    def set_stoich(self, stoichiometry: float):
        self.stoichiometry = stoichiometry

    def get_stoich(self):
        return self.stoichiometry

    def set_kmer(self, kmer: int):
        self.kmer = kmer

    def get_kmer(self):
        return self.kmer

    def set_ppi(self, ppi: str):
        if ppi not in self.ppis:
            self.ppis.append(ppi)

    def get_ppis(self):
        return self.ppis


class Complex:

    def __init__(self, ident: str, name: str):
        self.ident = ident
        self.name = name
        self.name_cd = None
        self.involved_molecules = []
        self.cat_reactions = []
        self.sclass = None
        self.stoichiometry = 1.
        self.definition = False
        self.structure = False
        self.position = False
        self.int_id = None
        self.ppis = []
        self.compartment = None

    def set_internal_id(self, int_id):
        self.int_id = int_id

    def get_internal_id(self):
        return self.int_id

    def get_element_by_name(self, name: str):
        for element in self.involved_molecules:
            if element.get_id() == name:
                return element
        return False

    def sort_molecules(self):
        self.involved_molecules = sorted(self.involved_molecules,
                                         key=lambda x: x.get_name())

    def get_id(self):
        return self.ident

    def get_name(self):
        return self.name

    def set_compartment(self, compartment):
        self.compartment = compartment

    def get_compartment(self):
        return self.compartment

    def add_molecule(self, Molecule: Molecule):
        self.involved_molecules.append(Molecule)

    def get_molecules(self):
        return self.involved_molecules

    def set_position(self, position: int):
        self.position = position

    def get_position(self):
        return self.position

    def set_class(self, sclass: str):
        self.sclass = sclass

    def get_class(self):
        return self.sclass

    def set_name_cd(self, name_cd: str):
        self.name_cd = name_cd

    def get_name_cd(self):
        return self.name_cd

    def set_cat_reactions(self, cat_reactions):
        self.cat_reactions.append(cat_reactions)

    def get_cat_reactions(self):
        return self.cat_reactions

    def set_stoich(self, stoichiometry: float):
        self.stoichiometry = stoichiometry

    def get_stoich(self):
        return self.stoichiometry

    def set_definition(self, ComplexDef: ComplexDef):
        self.definition = ComplexDef

    def get_definition(self):
        return self.definition

    def set_structure(self, structure: str):
        self.structure = structure

    def get_structure(self):
        return self.structure

    def set_ppi(self, ppi: str):
        if ppi not in self.ppis:
            self.ppis.append(ppi)

    def get_ppis(self):
        return self.ppis


class Rule:
    def __init__(self, ident: str, name: str):
        self.ident = ident
        self.name = name
        self.lhs = []
        self.lhsp = []
        self.rhs = []
        self.rxn_type = False
        self.mod_type = False
        self.reversible = False
        self.modifiers = []
        self.modifier_name = False
        self.deg = []
        self.synth = []
        self.modification_types = {}
        self.single_substrates = []
        self.single_products = []
        self.temp_ident = False
        self.mod_targets = []

    def set_temp_ident(self, temp_ident):
        self.temp_ident = temp_ident

    def get_temp_ident(self):
        return self.temp_ident

    def get_id(self):
        return self.ident

    def set_id(self, ident: str):
        self.ident = ident

    def get_name(self):
        return self.name

    def set_name(self, name: str):
        self.name = name

    def set_deg(self, deg_elem):
        if deg_elem not in self.deg:
            self.deg.append(deg_elem)

    def get_deg(self):
        return self.deg

    def add_single_substrate(self, substrate):
        self.single_substrates.append(substrate)

    def get_single_substrates(self):
        return self.single_substrates

    def add_single_product(self, product):
        self.single_products.append(product)

    def get_single_products(self):
        return self.single_products

    def sort_entities(self):
        self.single_substrates = sorted(self.single_substrates,
                                        key=lambda x: x.get_name())
        self.single_products = sorted(self.single_products,
                                      key=lambda x: x.get_name())

    def set_synth(self, synth_elem):
        if synth_elem not in self.synth:
            self.synth.append(synth_elem)

    def get_synth(self):
        return self.synth

    def set_left_hand_molecule(self, Molecule: Molecule):
        self.lhs.append(Molecule)

    def set_right_hand_molecule(self, Molecule: Molecule):
        self.rhs.append(Molecule)

    def get_left_hand_side(self):
        return self.lhs

    def get_right_hand_side(self):
        return self.rhs

    def set_direction(self, reversible: bool):
        self.reversible = reversible

    def get_direction(self):
        return self.reversible

    def set_modifier(self, Molecule: Molecule):
        self.modifiers.append(Molecule)

    def set_modifier_name(self, name: str):
        self.modifier_name = name

    def get_modifier_name(self):
        return self.modifier_name

    def get_modifiers(self):
        return self.modifiers

    def set_reaction_type(self, reaction_type: str):
        self.reaction_type = reaction_type

    def get_reaction_type(self):
        return self.reaction_type

    def set_modification_type(self, modifier, modification_type):
        self.modification_types[modifier] = modification_type

    def get_modification_types(self):
        return self.modification_types

        
