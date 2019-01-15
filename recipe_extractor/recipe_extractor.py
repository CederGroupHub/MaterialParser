# coding=utf-8
__author__ = "Olga Kononova"
__maintainer__ = "Olga Kononova"
__email__ = "0lgaGkononova@yandex.ru"
__version__ = "4.0"

import itertools
import re
from pprint import pprint
import sympy
from sympy.abc import _clash
from material_parser.material_parser import MaterialParser


# noinspection PyBroadException
class RecipeExtractor:
    def __init__(self, verbose=False, pubchem_lookup=False):
        print('RecipeExtractor version 4.5')
        self.__mp = MaterialParser(pubchem_lookup=pubchem_lookup, verbose=verbose)
        self.__verbose = verbose
        self.__pubchem = pubchem_lookup

        self.__greek_letters = [chr(i) for i in range(945, 970)]

        self.__list_of_elements_1 = ['H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U']
        self.__list_of_elements_2 = ['He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'Cr',
                                     'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',
                                     'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'Xe',
                                     'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
                                     'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
                                     'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
                                     'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
                                     'Fl', 'Lv']

    def get_composition(self, abstract_materials_, synthesis_materials, abstract=None, syn_paragraph=None):
        """
        main method to convert list of materials into composition
        based on output of MER
        :param abstract_materials_: <dict> targets: list of materials
                                            precursors: list of materials
                                            others: list of materials
        :param synthesis_materials: <dict> targets: list of materials
                                            precursors: list of materials
                                            others: list of materials
        :param abstract: <list> of <str> abstract sentences
        :param syn_paragraph: <list> of <str> paragraph sentences
        :return: output_structure: <dict> targets: list of material structures (output of mp.parse_material)
                                          precursors: list of material structures
                                          targets: list of material structures (output of mp.parse_material)
                                          others: list of material structures
                fails: <dict> targets: list of failed material strings
                              precursors: list of failed material strings
                              targets: list of failed material strings
                              others: list of failed material strings
                abbreviations: <dict> abbreviation: corresponding name

        """

        if abstract is None:
            abstract = []
        output_structure = {}

        fails = {'targets': [],
                 'precursors': [],
                 'abstracts': [],
                 'others': []}

        abstract_materials = list(set([m for m in abstract_materials_['targets']]))
        targets = list(set([m for m in synthesis_materials['targets']]))
        precursors = list(set([m for m in synthesis_materials['precursors']]))
        others = list(set([m for m in synthesis_materials['others']+abstract_materials_['others'] if m not in precursors + targets]))

        if self.__verbose:
            print('MER found materials')
            print('\tfrom abstract:', abstract_materials)
            print('\ttargets from synthesis paragaph', targets)
            print('\tprecursors from synthesis paragaph', precursors)
            print('\tother materials from synthesis paragaph and abstract', others)

        # Building abbreviations dictionary

        materials = list(set(abstract_materials + targets + others))
        abbreviations = self.__mp.build_abbreviations_dict(materials, abstract + syn_paragraph)

        if self.__verbose:
            print('Abbreviations dictionary:')
            pprint(abbreviations)

        def substitute_abbreviations(materials_list=None, abbreviations_dict=None):
            if abbreviations_dict is None:
                abbreviations_dict = {}
            if materials_list is None:
                materials_list = []
            materials_list_upd = []
            for _ in materials_list:
                if _ in abbreviations_dict:
                    materials_list_upd.append(abbreviations_dict[_])
                else:
                    materials_list_upd.append(_)

            return list(set(materials_list_upd))

        abstract_materials_upd_2 = substitute_abbreviations(abstract_materials, abbreviations)
        targets_upd_2 = substitute_abbreviations(targets, abbreviations)
        precursors_upd_2 = substitute_abbreviations(precursors, abbreviations)
        others_upd_2 = substitute_abbreviations(others, abbreviations)

        if self.__verbose:
            print('After abbreviations substitution:')
            print('\tabstract:', abstract_materials_upd_2)
            print('\ttargets:', targets_upd_2)
            print('\tprecursors:', precursors_upd_2)
            print('\tothers:', others_upd_2)

        # Resolving list of materials in name
        def get_list(material_name):
            if self.__mp.is_materials_list(material_name):
                m_list = self.__mp.reconstruct_list_of_materials(material_name)
                materials_list = []
                for mat, v in m_list:
                    f = self.__mp.reconstruct_formula(mat, v)
                    f = mat if f == '' else f
                    materials_list.append(f)
            else:
                materials_list = [material_name]

            return materials_list

        abstract_materials_upd_3 = list(set([m for material in abstract_materials_upd_2 for m in get_list(material)]))
        targets_upd_3 = list(set([m for material in targets_upd_2 for m in get_list(material)]))
        precursors_upd_3 = list(set([m for material in precursors_upd_2 for m in get_list(material)]))
        others_upd_3 = list(set([m for material in others_upd_2 for m in get_list(material)]))

        if self.__verbose:
            print('After reconstruction of materials list:')
            print('\tabstract:', abstract_materials_upd_3)
            print('\ttargets:', targets_upd_3)
            print('\tprecursors:', precursors_upd_3)
            print('\tothers:', others_upd_3)

        # Parsing chemical formulas of materials
        if self.__verbose:
            print('Extracting chemical composition:')

        def chemical_structure(materials_list=None):
            materials_structures = []

            for mat in materials_list:
                t_struct = self.__mp.parse_material(mat)
                if t_struct['composition'] != {}:
                    materials_structures.append(t_struct)

            return materials_structures

        abstract_materials_struct = chemical_structure(abstract_materials_upd_3)
        targets_struct = chemical_structure(targets_upd_3)
        precursors_struct = chemical_structure(precursors_upd_3)
        others_struct = chemical_structure(others_upd_3)

        # Resolving variables
        def resolve_variables(material_structure):

            sentences = abstract+syn_paragraph
            for _ in abstract_materials + targets + precursors + others:
                if any(s in _ for s in ['⩽', '<', '=']):
                    sentences = ['(' + _ + ')'] + sentences

            for var in material_structure['fraction_vars']:
                values = self.get_values(var, sentences)
                material_structure['fraction_vars'][var] = values

            for var in material_structure['elements_vars']:
                values = self.get_values(var, syn_paragraph, elements=True)
                if len(values) == 0:
                    values = self.get_values(var, abstract, elements=True)
                material_structure['elements_vars'][var] = values

            return material_structure

        abstract_materials_struct = [resolve_variables(m) for m in abstract_materials_struct]
        targets_struct = [resolve_variables(m) for m in targets_struct]
        others_struct = [resolve_variables(m) for m in others_struct]

        for material in precursors_struct:
            for x in material['fraction_vars'].keys():
                for target in targets_struct:
                    if x in target['fraction_vars']:
                        material['fraction_vars'][x] = target['fraction_vars'][x]

            for el in material['elements_vars'].keys():
                for target in targets_struct:
                    if el in target['elements_vars']:
                        material['elements_vars'][el] = target['elements_vars'][el]

        if self.__verbose:
            print('Materials structures from abstract:')
            pprint(abstract_materials_struct)
            print('Materials structures from targets:')
            pprint(targets_struct)
            print('Materials structures from precursors:')
            pprint(precursors_struct)
            print('Materials structures from other materials:')
            pprint(others_struct)

        # assigning possible abbreviations:
        for material in abstract_materials_struct:
            material['is_abbreviation_like'] = self.is_abbreviation_like(material)

        for material in targets_struct:
            material['is_abbreviation_like'] = self.is_abbreviation_like(material)

        for material in others_struct:
            material['is_abbreviation_like'] = self.is_abbreviation_like(material)

        if not targets_struct:
            targets_struct = abstract_materials_struct

        final_targets = [m for m in targets_struct if not m['is_abbreviation_like'] and m['composition'] != {}]

        output_structure['targets'] = final_targets
        output_structure['precursors'] = precursors_struct
        output_structure['abstract'] = abstract_materials_struct
        output_structure['others'] = others_struct

        return output_structure, fails, abbreviations

    ###############################################################################################################

    def substitute_elements(self, material_structure):
        """
        substituting values for elements variables into formula
        :param material_structure: <dict> output of mp.parse_material() with filled "element_vars"
        :return: list of structures derived from input with substitution of all element_vars
        """

        new_materials_array = []
        elements_array = self.__get_substitutions_array(material_structure['elements_vars'])

        for subs in elements_array:
            material_composition = {}
            for m in material_structure['composition'].keys():
                material_composition[m] = material_structure['composition'][m].copy()
                material_composition[m]['composition'] = material_structure['composition'][m]['composition'].copy()
                for var, val in subs.items():
                    material_composition[m]['composition'][val] = material_composition[m]['composition'][var]
                    del material_composition[m]['composition'][var]

            new_materials_array.append(dict(
                substitution=subs,
                updated_material=material_composition,
                fraction_vars=material_structure['fraction_vars'].copy()
            ))

        return new_materials_array

    def substitute_fraction(self, material_structure):
        """
        substituting values for elements fractions variables into formula
        :param material_structure: <dict> output of mp.parse_material() with filled "fraction_vars"
        :return: list of structures derived from input with substitution of all fraction_vars
        """

        new_materials_array = []
        fraction_variables = {x: v['values'] for x, v in material_structure['fraction_vars'].items()}
        fractions_array = self.__get_substitutions_array(fraction_variables)

        for subs in fractions_array:
            material_composition = {}
            for m in material_structure['composition'].keys():
                material_composition[m] = material_structure['composition'][m].copy()
                material_composition[m]['composition'] = material_structure['composition'][m]['composition'].copy()

                obtained = True
                for el, stoich in material_composition[m]['composition'].items():
                    for var, val in subs.items():
                        stoich = stoich.replace(var, str(val))
                    try:
                        stoich = round(float(eval(stoich)), 3)
                        if stoich < 0:
                            obtained = False
                    except:
                        obtained = False

                    if obtained:
                        material_composition[m]['composition'][el] = str(stoich)
                    else:
                        material_composition[m]['composition'][el] = material_structure[m]['composition'][el]

            new_materials_array.append(dict(
                updated_material=material_composition,
                substitution=subs,
            ))

        return new_materials_array

    def get_values(self, var, sents, elements=False):

        if not elements:
            values = dict(values=[], max_value=None, min_value=None)
            i = 0
            while not (len(values['values']) != 0 or values['max_value'] is not None) and i < len(sents):
                values = self.__mp.get_stoichiometric_values(var, sents[i])
                #print ('->', sents[i])
                #print ('->', values)
                i += 1
        else:
            values = []
            i = 0
            while len(values) == 0 and i < len(sents):
                values = self.__mp.get_elements_values(var, sents[i].strip('., '))
                i += 1

        return values

    def __get_substitutions_array(self, subs_dict):

        """
        Generate combinations of different variables values
        I.e. if 'x' = [0.1, 0.2] and 'y' = [0.5, 0.6], then outputs: [
        {'x': 0.1, 'y': 0.5}, {'x': 0.1, 'y': 0.6}, {'x': 0.2, 'y': 0.5},  {'x': 0.2, 'y': 0.6}]
        :return:
        """

        subs_array = []

        l_dict = len(subs_dict)
        t_array = [dict(var=k, val=v) for k, vs in subs_dict.items() for v in vs]

        for comb in itertools.combinations(range(0, len(t_array)), l_dict):
            s = ''.join([t_array[i]['var'] for i in comb])
            if len(s) == len(set(s)):
                t_dict = {}
                for i in comb:
                    t_dict[t_array[i]['var']] = t_array[i]['val']
                subs_array.append(t_dict)

        return subs_array

    def is_abbreviation_like(self, structure):
        if len(structure['elements_vars']) > 1 and \
                all(len(values) == 0 for values in structure['elements_vars'].values()):
            return True

        # if all(el[0].isupper() and len(el) == 1 and value in ['1', '1.0'] for el, value in
        #        structure['composition'].items()) \
        #         and structure['composition'] != {}:
        #     return True

        # if all(el[0].isupper() and len(el) == 1 and v in ['1', '1.0'] for m, compos in structure['composition'].items()
        #        for el, v in compos['composition'].items()):
        #     return True

        # if any(all(el[0].isupper() and len(el) == 1 and v in ['1', '1.0'] for el, v in compos['composition'].items())
        #        for compos in structure['mixture'].values()):
        #     return True

        if any(all(el[0].isupper() and len(el) == 1 and v in ['1', '1.0'] for el, v in compos['composition'].items())
               for m, compos in structure['composition'].items()):
            return True

        return False

    def reconstruct_formula(self, init_formula, data):
        new_formula = ''

        r = '[A-Z]{1}[a-wyz]{0,1}'
        elements = re.findall(r, init_formula)
        if len(elements) == len(data['compos']):
            if all(el in data['compos'] for el in elements):
                new_formula = ''.join([el + self.__cast_stoichiometry(data['compos'][el]) for el in elements])
            else:
                for el in elements:
                    if el not in self.__list_of_elements_2 and len(el) == 2 and el not in data['subs']:
                        el = el[0]
                    if el in data['compos']:
                        new_formula = new_formula + el + self.__cast_stoichiometry(data['compos'][el])
                    elif el in data['subs']:
                        el_ = data['subs'][el]
                        new_formula = new_formula + el_ + self.__cast_stoichiometry(data['compos'][el_])
                    else:
                        new_formula = ''

        new_formula = new_formula.replace('1.0', '')

        return new_formula

    def __is_precursor(self, p_composition, t_composition):

        for l in self.__greek_letters:
            _clash[l] = sympy.Symbol(l)

        if any(not sympy.sympify(v, _clash).is_Number for v in p_composition.values()):
            return False

        precursors__me = [p for p in p_composition.keys() if p not in ['0', 'H', 'N', 'C', 'P']]
        target__me = [p for p in t_composition.keys() if p not in ['0', 'H', 'N', 'C', 'P']]
        if len(precursors__me) > len(target__me):
            return False

        if any(el in p_composition for el in t_composition.keys() if el not in ['O', 'H', 'C']) \
                and t_composition.keys() != p_composition.keys():
            return True

        return False

    def __cast_stoichiometry(self, value):

        value = float(value)
        if value == 1.0:
            return ''
        if value * 1000 % 1000 == 0.0:
            return str(int(value))

        return str(value)

