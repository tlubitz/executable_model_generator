#!/usr/bin/env python
from molecule import Molecule, Rule


class Model:
    '''
    A Model that comprises all information extracted from an CD SBML file.
    '''
    def __init__(self, ident: str, name: str):
        self.ident = ident
        self.name = name
        self.parameters = []
        self.molecules = []
        self.initial_conditions = []
        self.complexes = []
        self.rules = []
        self.simple_molecules = []

    def get_ident(self):
        return self.ident

    def get_name(self):
        return self.name

    def add_parameter(self, Parameter):
        self.parameters.append(Parameter)

    def get_parameters(self):
        return self.parameters

    def add_molecule(self, Molecule: Molecule):
        self.molecules.append(Molecule)

    def get_molecules(self):
        return self.molecules

    def add_initial_condition(self, InitialCondition):
        self.initial_conditions.append(InitialCondition)

    def get_initial_conditions(self):
        return self.initial_conditions

    def add_complex(self, Complex):
        self.complexes.append(Complex)

    def get_complexes(self):
        return self.complexes

    def remove_complex(self, complex_id: str):
        try: self.complex.remove(complex_id)
        except: pass

    def add_rule(self, Rule: Rule):
        self.rules.append(Rule)

    def get_rules(self):
        return self.rules

    def remove_rule(self, ident: str):
        for rule in self.rules:
            if rule.get_id() == ident:
                self.rules.remove(rule)

    def add_simple_molecule(self, s_molecule: str):
        if s_molecule not in self.simple_molecules:
            self.simple_molecules.append(s_molecule)

    def get_simple_molecules(self):
        return self.simple_molecules


class Parameter:
    '''
    A parameter from the CD SBML file which has a name and a value
    '''
    def __init__(self, name: str, value: float):
        self.name = name
        self.value = value

    def get_name(self):
        return self.name

    def get_value(self):
        return self.value


class InitialCondition:
    '''
    The initial condition for a given molecule of the model
    '''
    def __init__(self, Molecule: Molecule, value: float):
        self.molecule = Molecule
        self.value = value

    def get_value(self):
        return self.value

    def get_molecule(self):
        return self.molecule
