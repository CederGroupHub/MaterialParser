from collections import OrderedDict
import material_parser.core.chemical_sets as cs
from material_parser.core.constants import FORMULA_SYMBOLS_SET, AMOUNTS_SYMBOLS_SET

from pprint import pprint


class Composition:
    def __init__(self, composition=OrderedDict()):
        self._composition = composition

    def __repr__(self):
        return "Elements({composition})".format(composition={e: v for e,v in self._composition.items()})

    def __str__(self):
        return {(e,v) for e,v in self._composition}

    def __getitem__(self, item):
        if item not in self._composition:
            raise ValueError('No item {}'.format(item))
        return self._composition[item]

    def __setitem__(self, key, value):
        if not isinstance(value, str):
            raise TypeError('Expected {} to be an str'.format(value))
        self._composition[key] = value

    @property
    def data(self):
        return self._composition

    @data.setter
    def data(self, value):
        if not isinstance(value, OrderedDict):
            raise TypeError('Expected {} to be an OrderedDict'.format(value))
        self._composition = OrderedDict([(k, v) for k, v in value.items()])


class Compound:
    def __init__(self, formula="", amount=1, elements=OrderedDict(), species=OrderedDict()):
        self._formula = formula
        self._amount = amount
        self._elements = Composition(elements)
        self._species = Composition(species)

    def __set__(self, obj, data):
        if not isinstance(data, dict):
            raise TypeError('Expected {} to be a dict'.format(data))
        obj.formula = data.get("formula", "")
        obj.amount = data.get("amount", "")
        obj.elements = data.get("elements", OrderedDict())
        obj.species = data.get("species", OrderedDict())

    @property
    def amount(self):
        return self._amount

    @amount.setter
    def amount(self, value):
        if not isinstance(value, str):
            raise TypeError('Expected {} to be an str'.format(value))
        if any(c not in AMOUNTS_SYMBOLS_SET for c in value):
            raise ValueError('Expected {amount} symbols to be one of {sym_set}'
                             .format(amount=value, sym_set=AMOUNTS_SYMBOLS_SET))
        self._amount = value

    @property
    def formula(self):
        return self._formula

    @formula.setter
    def formula(self, value):
        if not isinstance(value, str):
            raise TypeError('Expected {} to be an str'.format(value))
        if any(c not in FORMULA_SYMBOLS_SET for c in value):
            raise ValueError('Expected {amount} symbols to be one of {sym_set}'
                             .format(amount=value, sym_set=FORMULA_SYMBOLS_SET))
        self._formula = value

    @property
    def elements(self):
        return self._elements.data

    @elements.setter
    def elements(self, value):
        if not isinstance(value, dict):
            raise TypeError('Expected {} to be a dict'.format(value))
        self._elements.data = value

    @property
    def species(self):
        return self._species.data

    @species.setter
    def species(self, value):
        if not isinstance(value, dict):
            raise TypeError('Expected {} to be a dict'.format(value))
        self._species.data = value

    def __repr__(self):
        return "Compound({data_dict})".format(data_dict=self.to_dict())

    def __str__(self):
        return self.to_dict()

    def to_dict(self):
        return dict(formula=self._formula,
                    amount=self._amount,
                    elements=self._elements.data,
                    species=self._species.data)


# class Variables:
#     def __init__(self, variables={}):
#         self._variables = variables
#
#     def __set__(self, obj, data):
#         if not isinstance(data, dict):
#             raise TypeError('Expected {} to be a dict'.format(data))
#         obj.data = data
#
#     def __get__(self, instance, owner=None):
#         return instance.data
#
#     def __repr__(self):
#         return "Variables({})".format(self._variables)
#
#     def __str__(self):
#         return "({})".format(self._variables)
#
#     def __getitem__(self, item):
#         if item not in self._variables:
#             raise ValueError('No item {}'.format(item))
#         return self._variables[item]
#
#     def __setitem__(self, key, value):
#         if not isinstance(value, list):
#             raise TypeError('Expected {} to be a list'.format(value))
#         self._composition[key] = [v for v in value]
#
#     @property
#     def data(self):
#         return self._variables
#
#     @data.setter
#     def data(self, value):
#         if not isinstance(value, dict):
#             raise TypeError('Expected {} to be a dict'.format(value))
#         if any(not isinstance(k, str) for k in value.keys()):
#             raise TypeError('Expected keys to be an str.')
#         self._variables = {k: [v for v in vals] for k, vals in value.items()}


