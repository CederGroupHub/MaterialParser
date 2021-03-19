from collections import OrderedDict
import chemical_sets as cs
from data_descriptors import Variables, Compound
from constants import DEFICIENCY_CHARS, FORMULA_SYMBOLS_SET


class ChemicalStructure:
    def __init__(self, init_string):
        self._material_string = init_string
        self._material_name = ""
        self._material_formula = ""

        self._additives = []
        self._phase = ""
        self._oxygen_deficiency = ""
        self._is_acronym = False

        self._amounts_x = Variables()
        self._elements_x = Variables()

        self._composition = []

    def __set__(self, obj, data):
        if not isinstance(data, dict):
            raise TypeError('Expected {} to be a dict'.format(data))
        obj.material_name = data.get("material_name", "")
        obj.material_formula = data.get("material_formula", "")
        obj.additives = data.get("additives", [])
        obj.phase = data.get("phase", "")
        obj.oxygen_deficiency = data.get("oxygen_deficiency", "")
        obj.is_acronym = data.get("is_acronym", False)
        obj.amounts_x = data.get("amounts_x", {})
        obj.elements_x = data.get("elements_x", {})
        obj.composition = data.get("composition", [])

    @property
    def material_name(self):
        return self._material_name

    @material_name.setter
    def material_name(self, value):
        if not isinstance(value, str):
            raise TypeError('Expected {} to be an str'.format(value))
        self._material_name = value

    @property
    def material_formula(self):
        return self._material_formula

    @material_formula.setter
    def material_formula(self, value):
        if not isinstance(value, str):
            raise TypeError('Expected {} to be an str'.format(value))
        if any(c not in FORMULA_SYMBOLS_SET for c in value):
            # for c in value:
            #     if c not in FORMULA_SYMBOLS_SET:
            #         print ("--->", c)
            raise ValueError('Expected {value} symbols to be one of {sym_set}'
                             .format(value=value, sym_set=FORMULA_SYMBOLS_SET))
        self._material_formula = value

    @property
    def additives(self):
        return self._additives

    @additives.setter
    def additives(self, value):
        if not isinstance(value, list):
            raise TypeError('Expected {} to be a list'.format(value))
        if any(not isinstance(a, str) for a in value):
            raise TypeError('Expected {} to have items of an str'.format(value))
        self._additives = [a for a in value]

    @property
    def phase(self):
        return self._phase

    @phase.setter
    def phase(self, value):
        if not isinstance(value, str):
            raise TypeError('Expected {} to be an str'.format(value))
        self._phase = value

    @property
    def oxygen_deficiency(self):
        return self._oxygen_deficiency

    @oxygen_deficiency.setter
    def oxygen_deficiency(self, value):
        if not isinstance(value, str):
            raise TypeError('Expected {} to be an str'.format(value))
        if len(value) > 1:
            raise ValueError('Expected {oxy_def} to be a symbol'.format(value))
        if value not in DEFICIENCY_CHARS:
            raise ValueError('Expected {oxy_def} symbols to be one of {sym_set}'
                             .format(oxy_def=value, sym_set=DEFICIENCY_CHARS))
        self._oxygen_deficiency = value

    @property
    def is_acronym(self):
        return self._is_acronym

    @is_acronym.setter
    def is_acronym(self, value):
        if not isinstance(value, bool):
            raise TypeError('Expected {value} to be a boolean'.format(value))
        self._is_acronym = value

    @property
    def elements_x(self):
        return self._elements_x.data

    @elements_x.setter
    def elements_x(self, value):
        if not isinstance(value, dict):
            raise TypeError('Expected {} to be a dict'.format(value))
        if any(not isinstance(k, str) for k in value.keys()):
            raise TypeError('Expected {} to have keys od a str'.format(value))
        if any(not isinstance(v, list) for v in value.values()):
            raise TypeError('Expected {} to have values of alist'.format(value))
        self._elements_x.data = value

    @property
    def amounts_x(self):
        return self._amounts_x.data

    @amounts_x.setter
    def amounts_x(self, value):
        if not isinstance(value, dict):
            raise TypeError('Expected {} to be a dict'.format(value))
        if any(not isinstance(k, str) for k in value.keys()):
            raise TypeError('Expected {} to have keys of a str'.format(value))
        if any(not isinstance(v, list) for v in value.values()):
            raise TypeError('Expected {} to have values of a list'.format(value))
        self._amounts_x.data = value

    @property
    def composition(self):
        return [c for c in self._composition]

    @composition.setter
    def composition(self, value):
        if not isinstance(value, list):
            raise TypeError('Expected {} to be a list'.format(value))
        if any(not isinstance(c, dict) for c in value):
            raise TypeError('Expected {} to have items of type dict'.format(value))
        self._composition = [Compound(formula=v.get("formula", ""),
                                      amount=v.get("amount", ""),
                                      elements=v.get("elements", OrderedDict()),
                                      species=v.get("species", OrderedDict()))
                             for v in value]

    def element_structure(self, element):
        self._material_name = cs.element2name[element]
        self._material_formula = element
        self._composition = [Compound(element, 1, [(element, 1)], [(element, 1)])]
        return self

    def diatomic_molecule_structure(self, element):
        self._material_name = cs.element2name[element]
        self._material_formula = element
        self._composition = [Compound(element, 1, [(element, 2)], [(element, 1)])]
        return self

    def water_structure(self):
        self._material_name = "water"
        self._material_formula = "H2O"
        self._composition = [Compound("H2O", 1, [("H", 2), ("O", 1)], [("H2O", 1)])]
        return self

    def to_dict(self):
        return dict(material_string=self._material_string,
                    material_name=self._material_name,
                    material_formula=self._material_formula,
                    additives=self._additives,
                    phase=self._phase,
                    oxygen_deficiency=self._oxygen_deficiency,
                    is_acronym=self._is_acronym,
                    amounts_x=self._amounts_x.data,
                    elements_x=self._elements_x.data,
                    composition=[c.to_dict() for c in self._composition])
