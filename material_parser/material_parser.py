# coding=utf-8

__author__ = "Olga Kononova"
__maintainer__ = "Olga Kononova"
__email__ = "0lgaGkononova@yandex.ru"
__version__ = "3.0"

import regex as re
import collections
import sympy as smp
from sympy.abc import _clash
import pubchempy as pcp
import os
import json


# noinspection PyBroadException
class MaterialParser:
    def __init__(self, verbose=False, pubchem_lookup=False, fails_log=False):
        print('MaterialParser version 3.6')

        self.__list_of_elements_1 = ['H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U']
        self.__list_of_elements_2 = ['He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'Cr',
                                     'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',
                                     'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'Xe',
                                     'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
                                     'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
                                     'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
                                     'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
                                     'Fl', 'Lv']
        self.__greek_letters = [chr(i) for i in range(945, 970)]

        self.__filename = os.path.dirname(os.path.realpath(__file__))

        self.__ions = json.loads(open(os.path.join(self.__filename, 'rsc/ions_dictionary.json')).read())
        self.__anions = {ion['c_name']: {'valency': ion['valency'], 'e_name': ion['e_name'], 'n_atoms': ion['n_atoms']}
                         for ion in self.__ions['anions']}
        self.__cations = {ion['c_name']: {'valency': ion['valency'], 'e_name': ion['e_name'], 'n_atoms': ion['n_atoms']}
                          for ion in self.__ions['cations']}
        self.__chemicals = self.__ions['chemicals']
        self.__element2name = self.__ions['elements']

        self.__prefixes2num = {'': 1, 'mono': 1, 'di': 2, 'tri': 3, 'tetra': 4, 'pent': 5, 'penta': 5, 'hexa': 6,
                               'hepta': 7, 'octa': 8, 'nano': 9, 'ennea': 9, 'nona': 9, 'deca': 10, 'undeca': 11,
                               'dodeca': 12}
        self.__neg_prefixes = ['an', 'de', 'non']

        self.__rome2num = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10}

        self.__pubchem_dictionary = json.loads(open(os.path.join(self.__filename, 'rsc/pubchem_dict.json')).read())

        self.__fails_log = fails_log
        if fails_log:
            self.__pubchem_file = open(os.path.join(self.__filename, 'pubchem_log'), 'w')
            self.__pubchem_file.close()

        self.__pubchem = pubchem_lookup
        if pubchem_lookup:
            print ('Pubchem lookup is on! Will search for unknown materials name in PubChem DB.')

        self.__verbose = verbose

    ###################################################################################################################
    # Parsing material name
    ###################################################################################################################

    def parse_material(self, material_string_):
        """
        Main method to parse material string into chemical structure and convert chemical name into chemical formula
        :param material_string_: < str> material name/formula
        :return: dict(material_string: <str> initial material string,
                     material_name: <str> chemical name of material found in the string
                     material_formula: <str> chemical formula of material
                     dopants: <list> list of dopped materials/elements appeared in material string
                     phase: <str> material phase appeared in material string
                     hydrate: <float> if material is hydrate fraction of H2O
                     is_mixture: <bool> material is mixture/composite/alloy/solid solution
                     is_abbreviation: <bool> material is similar to abbreviation
                     fraction_vars: <dict> elements fraction variables: <list> values
                     elements_vars: <dict> elements variables: <list> values
                     composition: <dict> compound constitute of the material: its composition (element: fraction) and
                                                                            fraction of compound)
        """

        material_string = self.cleanup_name(material_string_)

        if self.__verbose:
            print ('After cleaning up string:')
            print (material_string_, '-->', material_string)

        dopants, material_string = self.get_dopants(material_string)
        if self.__verbose:
            print ('After dopants extraction:')
            print (material_string, 'with', dopants)

        material_string = material_string.lstrip(') -').rstrip('( ,.:;-±/δ+')
        material_name, material_formula, material_structure = self.split_material_name(material_string)

        if self.__verbose:
            print ('After material name parsing:')
            print (material_string, '-->', material_name, 'and', material_formula)

        # if material_string contains chemical formula
        if material_structure['composition'] != {}:
            output_structure = dict(
                material_string=material_string_,
                material_name=material_name,
                material_formula=material_formula,
                dopants=dopants,
                phase=material_structure['phase'],
                hydrate=material_structure['hydrate'],
                is_mixture=False,
                is_abbreviation_like=False,
                fraction_vars=material_structure['fraction_vars'],
                elements_vars=material_structure['elements_vars'],
                composition={
                    material_structure['formula']: dict(
                        composition=material_structure['composition'],
                        fraction='1'
                    )
                }
            )
            # noinspection PyBroadException
            try:
                output_structure['hydrate'] = float(output_structure['hydrate'])
            except:
                output_structure['hydrate'] = None

            return output_structure
        else:
            material_formula = ''

        # if material_string is chemical name reconstructing its formula
        if len(material_name) > 0:
            material_formula = self.reconstruct_formula(material_name)

            if material_formula == '':
                for m in [material_name,
                          material_name.lower(),
                          material_name.replace('-', ' '),
                          material_name.replace('-', ' ').lower()]:
                    if m in self.__pubchem_dictionary and material_formula == '':
                        material_formula = self.__pubchem_dictionary[m]

        if self.__verbose:
            print ('After looking up in dictionary:')
            print (material_string, '-->', material_name, 'and', material_formula)

        material_formula = material_string if material_formula == '' else material_formula
        # noinspection PyBroadException
        try:
            material_parts = self.split_material(material_formula)
        except:
            material_parts = [material_formula]

        if self.__verbose:
            print ('After splitting:')
            print (material_string, '-->', material_name, material_parts)

        output_structure = dict(
            material_string=material_string_,
            material_name=material_name,
            material_formula=material_formula.replace(' ', ''),
            phase='',
            hydrate='',
            dopants=dopants,
            is_mixture=len(material_parts) > 1,
            is_abbreviation_like=False,
            fraction_vars={},
            elements_vars={},
            composition={}
        )

        for m, f in material_parts:
            try:
                m = self.__check_parentheses(m)
                structure = self.get_structure_by_formula(m)
                output_structure['phase'] = structure['phase']
                output_structure['hydrate'] = structure['hydrate']
                output_structure['fraction_vars'].update(structure['fraction_vars'])
                output_structure['elements_vars'].update(structure['elements_vars'])
                if m == 'H2O':
                    output_structure['hydrate'] = f
                elif structure['composition'] != {}:
                    output_structure['composition'][structure['formula']] = dict(
                        composition=structure['composition'],
                        fraction=f
                    )
            except:
                output_structure['composition'][m] = dict(
                    composition={},
                    fraction=f
                )

        # converting hydrate fraction to float
        try:
            output_structure['hydrate'] = float(output_structure['hydrate'])
        except:
            output_structure['hydrate'] = None

        # substituting dopant into composition if it makes fractions to sum-up to integer
        dopant = dopants[0].strip(' ') if len(dopants) == 1 else ''
        #print('-->', dopant)
        if dopant != '':
            output_structure['material_formula'], output_structure['composition'] = \
                self.__substitute_dopant(dopant, material_formula, output_structure['composition'])

        # checking ions
        ions_set = set(self.__list_of_elements_2+self.__list_of_elements_1)-set(['H', 'N', 'O', 'Ar'])
        if material_formula.rstrip('0987654321+') in ions_set:
            output_structure['composition'] = {}


        if output_structure['composition'] == {} and self.__fails_log:
            with open(os.path.join(self.__filename, 'fails_log'), 'a') as f_pubchem:
                f_pubchem.write(material_name + '\n')

        return output_structure

    def get_structure_by_formula(self, formula):
        '''
        Parsing chemical formula in composition
        :param formula: <str> chemical formula
        :return: dict(formula: <str> formula string corresponding to obtained composition
                     composition: <dict> element: fraction
                     fraction_vars: <dict> elements fraction variables: <list> values
                     elements_vars: <dict> elements variables: <list> values
                     hydrate: <str> if material is hydrate fraction of H2O
                     phase: <str> material phase appeared in formula
                    )
        '''

        formula = formula.replace(' ', '')
        formula = formula.replace('−', '-')
        formula = formula.strip(' ')

        # is there any phase specified
        phase = ''
        if formula[0].islower():
            for m in re.finditer('([' + ''.join(self.__greek_letters) + ']*)-{0,1}(.*)', formula):
                if m.group(2) != '':
                    phase = m.group(1)
                    formula = m.group(2)

        # checking for hydrate
        hyd_i = formula.find('H2O') - 1
        hydrate_num = []
        while hyd_i > 0 and formula[hyd_i] not in ['·', '•', '-', '×', '⋅']:
            hydrate_num.append(formula[hyd_i])
            hyd_i -= 1
        hydrate = ''.join([c for c in reversed(hydrate_num)])
        if hyd_i > 0:
            formula = formula[:hyd_i]

        elements_variables = collections.defaultdict(str)
        stoichiometry_variables = collections.defaultdict(str)

        # convert fractions a(b-/+x) into a*b-/+a*x
        for m in re.findall('([0-9\.]+)(\([0-9\.a-z]+[-+]+[0-9\.a-z]+\))', formula):
            expr = str(smp.simplify(m[0] + '*' + m[1]))
            if expr[0] == '-':
                s_expr = re.split('\+', expr)
                expr = s_expr[1] + s_expr[0]
            expr = expr.replace(' ', '')
            formula = formula.replace(m[0] + m[1], expr, 1)

        # check for any weird syntax
        r = "\(([^\(\)]+)\)\s*([-*\.\da-z\+/]*)"
        for m in re.finditer(r, formula):
            if ',' in m.group(1):
                elements_variables['M'] = re.split(',', m.group(1))
                formula = formula.replace('(' + m.group(1) + ')' + m.group(2), 'M' + m.group(2), 1)
            if not m.group(1).isupper() and m.group(2) == '':
                formula = formula.replace('(' + m.group(1) + ')', m.group(1), 1)

        composition = self.__parse_formula(formula)

        # looking for variables in elements and stoichiometry
        for el, amt in composition.items():
            if el not in self.__list_of_elements_1 + \
                    self.__list_of_elements_2 + \
                    list(elements_variables.keys()) + ['□']:
                elements_variables[el] = []
            for var in re.findall('[a-z' + ''.join(self.__greek_letters) + ']', amt):
                stoichiometry_variables[var] = {}

        if 'R' in elements_variables and 'E' in elements_variables:
            elements_variables['RE'] = []
            del elements_variables['E']
            del elements_variables['R']
            composition['RE'] = composition['E']
            del composition['R']
            del composition['E']

        if 'A' in elements_variables and 'E' in elements_variables:
            elements_variables['AE'] = []
            del elements_variables['E']
            del elements_variables['A']
            composition['AE'] = composition['E']
            del composition['R']
            del composition['A']

        if 'M' in elements_variables and 'e' in stoichiometry_variables:
            elements_variables['Me'] = []
            del elements_variables['M']
            del stoichiometry_variables['e']
            c = composition['M'][1:]
            composition['Me'] = c if c != '' else '1.0'
            del composition['M']

        formula_structure = {'composition': composition,
                             'fraction_vars': stoichiometry_variables,
                             'elements_vars': elements_variables,
                             'hydrate': hydrate,
                             'phase': phase,
                             'formula': formula}

        return formula_structure

    def __parse_formula(self, init_formula):

        formula_dict = collections.defaultdict(str)

        formula_dict = self.__parse_parentheses(init_formula, "1", formula_dict)

        """
        refinement of non-variable values
        """
        incorrect = []
        for el, amt in formula_dict.items():
            formula_dict[el] = self.__simplify(amt)
            if any(len(c) > 1 for c in re.findall('[A-Za-z]+', formula_dict[el])):
                incorrect.append(el)

        for el in incorrect:
            del formula_dict[el]

        return formula_dict

    def __parse_parentheses(self, init_formula, init_factor, curr_dict):
        r = "\(((?>[^\(\)]+|(?R))*)\)\s*([-*\.\da-z\+/]*)"

        for m in re.finditer(r, init_formula):
            factor = "1"
            if m.group(2) != "":
                factor = m.group(2)

            factor = self.__simplify('(' + str(init_factor) + ')*(' + str(factor) + ')')

            unit_sym_dict = self.__parse_parentheses(m.group(1), factor, curr_dict)

            init_formula = init_formula.replace(m.group(0), '')

        unit_sym_dict = self.__get_sym_dict(init_formula, init_factor)
        for el, amt in unit_sym_dict.items():
            if len(curr_dict[el]) == 0:
                curr_dict[el] = amt
            else:
                curr_dict[el] = '(' + str(curr_dict[el]) + ')' + '+' + '(' + str(amt) + ')'

        return curr_dict

    def __get_sym_dict(self, f, factor):
        sym_dict = collections.defaultdict(str)
        r = "([A-Z□]{1}[a-z]{0,1})\s*([-*\.\da-z" + ''.join(self.__greek_letters) + "\+/]*)"

        el = ""
        amt = ""
        for m in re.finditer(r, f):
            """
            checking for correct elements names
            """
            el_bin = "{0}{1}".format(str(int(m.group(1)[0] in self.__list_of_elements_1 + ['M', '□'])), str(
                int(m.group(1) in self.__list_of_elements_1 + self.__list_of_elements_2 + ['Ln', 'M', '□'])))
            if el_bin in ['01', '11']:
                el = m.group(1)
                amt = m.group(2)
            if el_bin in ['10', '00']:
                el = m.group(1)[0]
                amt = m.group(1)[1:] + m.group(2)

            if len(sym_dict[el]) == 0:
                sym_dict[el] = "0"
            if amt.strip() == "":
                amt = "1"
            sym_dict[el] = '(' + sym_dict[el] + ')' + '+' + '(' + amt + ')' + '*' + '(' + str(factor) + ')'
            f = f.replace(m.group(), "", 1)
        if f.strip():
            return collections.defaultdict(str)
            # print("{} is an invalid formula!".format(f))

        """
        refinement of non-variable values
        """
        for el, amt in sym_dict.items():
            sym_dict[el] = self.__simplify(amt)

        return sym_dict

    ###################################################################################################################
    # Reconstruct formula / dictionary lookup
    ###################################################################################################################

    def split_material_name(self, material_string):
        '''
        Splitting material string into chemical name + chemical formula
        :param material_string: in form of "chemical name chemical formula"/"chemical name [chemical formula]"
        :return: name: <str> chemical name found in material string
                formula: <str> chemical formula found in material string
                structure: <dict> output of get_structure_by_formula()
        '''
        formula = ''
        structure = self.__empty_structure().copy()

        split = re.split('\s', material_string)

        if len(split) == 1:
            return material_string, formula, structure

        formula_e = split[-1].strip('[]')
        if re.match('(\s*\([IV,]+\))', formula_e):
            formula_e = ''
        try:
            structure_e = self.get_structure_by_formula(formula_e)
            composition_e = structure_e['composition']
        except:
            composition_e = {}
            formula_e = ''
            structure_e = self.__empty_structure().copy()

        formula_b = split[0].strip('[]')
        if re.match('(\s*\([IV,]+\))', formula_b):
            formula_b = ''
        try:
            structure_b = self.get_structure_by_formula(formula_b)
            composition_b = structure_b['composition']
        except:
            composition_b = {}
            formula_b = ''
            structure_b = self.__empty_structure().copy()

        if composition_e != {} and formula_e not in self.__list_of_elements_1 + self.__list_of_elements_2:
            split = split[:-1]
            structure = structure_e
            formula = formula_e
        elif composition_b != {} and formula_b not in self.__list_of_elements_1 + self.__list_of_elements_2:
            split = split[1:]
            structure = structure_b
            formula = formula_b
        else:
            formula = ''
            structure = self.__empty_structure().copy()

        name_terms = [p for p in split if
                      p.lower().strip('., -;:').rstrip('s') in self.__chemicals or 'hydrate' in p]

        if len(name_terms) > 0:
            name = ''.join([t + ' ' for t in split]).strip(' ')
        else:
            name = ''
            formula = ''
            structure = self.__empty_structure().copy()

        return name, formula, structure

    def reconstruct_formula(self, material_name, valency=''):
        """
        reconstructing chemical formula for simple chemical names anion + cation
        :param material_name: <str> chemical name
        :param valency: <str> anion valency
        :return: <str> chemical formula
        """

        output_formula = ''

        terms_list = []
        valency_list = []
        hydrate = ''
        cation_prefix_num = 0
        cation_data = {"c_name": "", "valency": [], "e_name": "", "n_atoms": 0}

        for t in material_name.split(' '):
            if t.strip('()') in self.__rome2num:
                valency_list.append(self.__rome2num[t.strip('()')])
                continue
            if 'hydrate' in t.lower():
                hydrate = t
                continue

            terms_list.append(t)

        t = ''.join([t + ' ' for t in terms_list]).lower().strip(' ')
        if t in self.__anions:
            return self.__anions[t]['e_name']
        if t in self.__cations:
            return self.__cations[t]['e_name']

        if len(terms_list) < 2:
            return output_formula

        if len(valency_list) > 1:
            print ('WARNING! Found many valencies per chemical name ' + material_name)
            print (valency_list)

        anion = terms_list.pop().lower().rstrip('s')

        if valency == '':
            valency_num = max(valency_list + [0])
        else:
            valency_num = self.__rome2num[valency.strip('()')]

        next_term = terms_list.pop()
        if 'hydrogen' in next_term.lower() and len(terms_list) != 0:
            anion = next_term + ' ' + anion
        else:
            terms_list += [next_term]

        _, anion_prefix_num, anion = self.__get_prefix(anion)

        if anion in self.__anions:
            anion_data = self.__anions[anion]
        else:
            return output_formula

        if len(terms_list) >= 2:
            return output_formula

        if len(terms_list) == 1:
            _, cation_prefix_num, cation = self.__get_prefix(terms_list[0])
            if cation.lower() in self.__cations:
                cation = cation.lower()
                cation_data = self.__cations[cation]
            elif cation in self.__element2name:
                cation_data = self.__cations[self.__element2name[cation]]
            else:
                return output_formula

        if len(cation_data['valency']) > 1 and valency_num != 0:
            if valency_num not in cation_data['valency']:
                print ('WARNING! Not common valency value for ' + material_name)
                print(cation_data['valency'])
                print(valency_num)
            cation_data['valency'] = [valency_num]

        output_formula = self.__build_formula(anion=anion_data,
                                              cation=cation_data,
                                              cation_prefix_num=cation_prefix_num,
                                              anion_prefix_num=anion_prefix_num)

        if hydrate != '':
            _, hydrate_prefix_num, hydrate = self.__get_prefix(hydrate)
            output_formula = output_formula + '·' + str(hydrate_prefix_num) + 'H2O'

        return output_formula

    def __build_formula(self, cation, anion, cation_prefix_num=0, anion_prefix_num=0):

        cation_stoich = 0
        anion_stoich = 0

        if anion_prefix_num + cation_prefix_num == 0 or anion_prefix_num * cation_prefix_num != 0:
            v_cation = abs(cation['valency'][0])
            v_anion = abs(anion['valency'][0])
            cm = self.__lcm(v_cation, v_anion)
            cation_stoich = cm // v_cation
            anion_stoich = cm // v_anion

        if anion_prefix_num != 0:
            anion_stoich = anion_prefix_num
            cation_stoich = anion_prefix_num * abs(anion['valency'][0]) // abs(cation['valency'][0])

        if cation_prefix_num != 0:
            cation_stoich = cation_prefix_num
            anion_stoich = cation_prefix_num * abs(cation['valency'][0]) // abs(anion['valency'][0])

        anion_name_el = anion['e_name']
        if anion_stoich > 1 and anion['n_atoms'] > 1:
            anion_name_el = '(' + anion_name_el + ')'

        cation_name_el = cation['e_name']
        if cation_stoich > 1 and cation['n_atoms'] > 1:
            cation_name_el = '(' + cation_name_el + ')'

        return "{0}{1}{2}{3}".format(cation_name_el, self.__cast_stoichiometry(cation_stoich), anion_name_el,
                                     self.__cast_stoichiometry(anion_stoich))

    def __get_prefix(self, material_name):

        pref = ''
        pref_num = 0
        material_name_upd = material_name

        for p in self.__prefixes2num.keys():
            if material_name.lower().find(p) == 0 and p != '':
                pref = p
                pref_num = self.__prefixes2num[p]
                material_name_upd = material_name_upd[len(p):].strip('-')

                if material_name_upd == 'xide':
                    material_name_upd = 'oxide'

        return pref, pref_num, material_name_upd

    ###################################################################################################################
    # Splitting mixtures
    ###################################################################################################################

    def split_material(self, material_name):
        """
        splitting mixture/composite/solid solution/alloy into compounds with fractions
        :param material_name: <str> material formula
        :return: <list> of <tuples>: (compound, fraction)
        """

        split = self.__split_name(material_name)
        l = 0
        while len(split) != l:
            l = len(split)
            split = [p for s in split for p in self.__split_name(s[0], s[1])]

        output = []
        for m, f in split:
            try:
                f = smp.simplify(f)
                if f.is_Number:
                    f = round(float(f), 3)
                f = str(f)
            except:
                f = '1'

            output.append((m, f))

        return output

    def __split_name(self, material_name_, init_fraction='1'):

        re_str = "(?<=[0-9\)])[-·∙\∗](?=[\(0-9])|(?<=[A-Z])[-·∙\∗](?=[\(0-9])|(?<=[A-Z\)])[-·∙\∗](?=[A-Z])|(?<=[" \
                 "0-9\)])[-·∙\∗](?=[A-Z]) "
        material_name = material_name_.replace(' ', '')

        if '(1-x)' == material_name[0:5]:
            material_name = material_name.replace('(x)', 'x')
            parts = re.findall('\(1-x\)(.*)[-+·∙\∗]x(.*)', material_name)
            parts = parts[0] if parts != [] else (material_name[5:], '')
            return [(parts[0].lstrip(' ·*'), '1-x'), (parts[1].lstrip(' ·*'), 'x')]

        parts = re.split(re_str, material_name)

        if len(parts) > 1:
            parts_upd = [p for part in parts for p in
                         re.split('(?<=[A-Z\)])[-·∙\∗](?=[xyz])|(?<=O[0-9\)]+)[-·∙\∗](?=[xyz])', part)]
        else:
            parts_upd = parts

        merged_parts = [parts_upd[0]]
        for m in parts_upd[1:]:
            if re.findall('[A-Z]', m) == ['O']:
                to_merge = merged_parts.pop() + '-' + m
                merged_parts.append(to_merge)
            else:
                merged_parts.append(m)

        composition = []
        for m in merged_parts:
            fraction = ''
            i = 0
            while i < len(m) and not m[i].isupper():
                fraction = fraction + m[i]
                i += 1
            fraction = fraction.strip('()')
            if fraction == '':
                fraction = '1'
            else:
                m = m[i:]

            fraction = '(' + fraction + ')*(' + init_fraction + ')'

            if m != '':
                composition.append((m, fraction))

        return composition

    def get_dopants(self, material_name):
        """
        resolving doped part in material string
        :param material_name: <str> material string
        :return: <list> of dopants,
                <str> new material name
        """
        new_material_name = material_name
        dopants = []

        # additions = ['doped', 'stabilized', 'activated','coated', 'modified']

        new_material_name = new_material_name.replace('codoped', 'doped')

        # checking for "doped with"
        for r in ['coated', 'activated', 'modified', 'stabilized', 'doped']:
            parts = [w for w in re.split(r + ' with', new_material_name) if w != '']
            if len(parts) > 1:
                new_material_name = parts[0].strip(' -+')
                dopants.append(parts[1].strip())

        # checking for element-doped prefix
        dopant = ''
        for r in ['coated', 'activated', 'modified', 'stabilized', 'doped']:
            parts = [w for w in re.split('(.*)[-\s]{1}' + r + ' (.*)', new_material_name) if w != '']
            if len(parts) > 1:
                new_material_name = parts.pop()
                dopants.extend(p for p in parts)

        if dopant != '':
            dopants.append(dopant)

        if '%' in new_material_name:
            new_material_name = new_material_name.replace('.%', '%')
            parts = re.split('[\-+:·]\s*[0-9\.]*\s*[vmolwtx\s]*\%', new_material_name)

        if len(parts) > 1:
            new_material_name = parts[0].strip(' -+')
            dopants.extend(d.strip(' ') for d in parts[1:] if d != '')

        for part in new_material_name.split(':'):
            if all(e.strip('x,+0987654321. ') in self.__list_of_elements_1 + self.__list_of_elements_2
                   for e in part.split(' ') if e != ''):
                dopants.append(part.strip(' '))
            else:
                new_material_name = part.strip(' ')

        return dopants, new_material_name

    ###################################################################################################################
    # Resolving abbreviations
    ###################################################################################################################

    def __is_abbreviation(self, word):
        if all(c.isupper() for c in re.sub('[0-9x\-\(\)\.]', '', word)) and len(re.findall('[A-NP-Z]', word)) > 1:
            return True

        return False

    def build_abbreviations_dict(self, materials_list, paragraph):
        """
        constructing dictionary of abbreviations appeared in material list
        :param paragraph: <list> of <str> list of sentences to look for abbreviations names
        :param materials_list: <list> of <str> list of materials entities
        :return: <dict> abbreviation: corresponding string
        """

        abbreviations_dict = {t: '' for t in materials_list if self.__is_abbreviation(t.replace(' ', '')) and t != ''}
        not_abbreviations = list(set(materials_list) - set(abbreviations_dict.keys()))

        # first find abreviations in current materials list
        for abbr in abbreviations_dict.keys():

            for material_name in not_abbreviations:
                if sorted(re.findall('[A-NP-Z]', abbr)) == sorted(re.findall('[A-NP-Z]', material_name)):
                    abbreviations_dict[abbr] = material_name

        # for all other abbreviations going through the paper text
        for abbr, name in abbreviations_dict.items():
            sents = ' '.join([s + ' ' for s in paragraph if abbr in s]).strip(' ').split(abbr)
            i = 0
            while abbreviations_dict[abbr] == '' and i < len(sents):
                sent = sents[i]
                for tok in sent.split(' '):
                    if sorted(re.findall('[A-NP-Z]', tok)) == sorted(re.findall('[A-NP-Z]', abbr)):
                        abbreviations_dict[abbr] = tok
                i += 1

        for abbr in abbreviations_dict.keys():
            parts = re.split('-', abbr)
            if all(p in abbreviations_dict for p in parts) and abbreviations_dict[abbr] == '' and len(parts) > 1:
                name = ''.join('(' + abbreviations_dict[p] + ')' + '-' for p in parts).rstrip('-')
                abbreviations_dict[abbr] = name

        empty_list = [abbr for abbr, name in abbreviations_dict.items() if name == '']
        for abbr in empty_list:
            del abbreviations_dict[abbr]

        return abbreviations_dict

    ###################################################################################################################
    # Methods to substitute variables
    ###################################################################################################################

    def __get_values(self, string, mode):
        values = []
        max_value = None
        min_value = None

        if len(string) == 0:
            return dict(values=[], max_value=None, min_value=None)

        # given range
        if mode == 'range':
            min_value, max_value = string[0]
            max_value = max_value.rstrip('., ')
            min_value = min_value.rstrip('., ')
            max_value = re.sub('[a-z]*', '', max_value)
            min_value = re.sub('[a-z]*', '', min_value)
            try:
                max_value = round(float(smp.simplify(max_value)), 4)
                min_value = round(float(smp.simplify(min_value)), 4) if min_value != '' else None
            except Exception as ex:
                max_value = None
                min_value = None
                template = "An exception of type {0} occurred when use sympy. Arguments:\n{1!r}."
                message = template.format(type(ex).__name__, ex.args)
                print(message)
            values = []

            return dict(values=values, max_value=max_value, min_value=min_value)

        # given list
        if mode == 'values':
            # values = re.split('[,\s]', string[0])
            values = re.split('[,\s]', re.sub('[a-z]+', '', string[0]))
            try:
                values = [round(float(smp.simplify(c.rstrip('., '))), 4) for c in values if
                          c.rstrip('., ') not in ['', 'and']]
                max_value = max(values) if values != [] else None
                min_value = min(values) if len(values) > 1 else None
            except Exception as ex:
                values = []
                max_value = None
                min_value = None
                template = "An exception of type {0} occurred when use sympy. Arguments:\n{1!r}"
                message = template.format(type(ex).__name__, ex.args)
                print(message)

        return dict(values=values, max_value=max_value, min_value=min_value)

    def get_stoichiometric_values(self, var, sentence):
        """
        find numeric values of var in sentence
        :param var: <str> variable name
        :param sentence: <str> sentence to look for
        :return: <dict>: max_value: upper limit
                        min_value: lower limit
                        values: <list> of <float> numeric values
        """
        values = dict(values=[], max_value=None, min_value=None)

        regs = [(var + '\s*=\s*([-]{0,1}[0-9\.\,/and\s]+)[\s\)\]\,]', 'values'),
                (var + '\s*=\s*([0-9\.]+)\s*[-–]\s*([0-9\.\s]+)[\s\)\]\,m\%]', 'range'),
                ('([0-9\.\s]*)\s*[<≤⩽]{0,1}\s*' + var + '\s*[<≤⩽>]{1}\s*([0-9\.\s]+)[\s\)\]\.\,]', 'range'),
                (var + '[a-z\s]*from\s([0-9\./]+)\sto\s([0-9\./]+)', 'range'),
                ]

        for r, m in regs:
            #print ('-->', r)
            if values['values'] == [] and values['max_value'] is None:
                r_res = re.findall(r, sentence.replace(' - ', '-'))
                #print ('-->', r_res)
                values = self.__get_values(r_res, m)
                #print ('-->', values)
            #print ('----')

        return values

    def get_elements_values(self, var, sentence):
        """
        find elements values for var in the sentence
        :param var: <str> variable name
        :param sentence: <str> sentence to look for
        :return: <list> of <str> found values
        """
        values = re.findall(var + '\s*[=:]{1}\s*([A-Za-z,\s]+)', sentence)
        if len(values) > 0:
            values = [c for c in re.split('[,\s]', values[0]) if
                      c in self.__list_of_elements_1 + self.__list_of_elements_2]

        return values

    ###################################################################################################################
    # Splitting list of materials
    ###################################################################################################################

    def is_materials_list(self, material_string):
        if any(a + 's' in material_string.lower() for a in self.__anions.keys()) and \
                any(w in material_string for w in ['and ', ',', ' of ']):
            return True

        return False

    def reconstruct_list_of_materials(self, material_string):
        """
        split material string into list of compounds when it's given in form cation + several anions
        :param material_string: <str>
        :return: <list> of <str> chemical names
        """

        parts = [p for p in re.split('[\s\,]', material_string) if p != '']

        anion = [(i, p[:-1]) for i, p in enumerate(parts) if p[:-1].lower() in self.__anions.keys()]
        cation = [(i, p) for i, p in enumerate(parts) if p.lower() in self.__cations.keys()
                  or p in self.__list_of_elements_1 + self.__list_of_elements_2]
        valencies = [(i - 1, p.strip('()')) for i, p in enumerate(parts) if p.strip('()') in self.__rome2num and i != 0]

        result = []
        if len(anion) == 1:
            for c_i, c in cation:

                if c in self.__element2name:
                    name = [self.__element2name[c]]
                else:
                    name = [c.lower()]
                valency = ''.join([v for v_i, v in valencies if v_i == c_i])
                if valency != '':
                    name.append('(' + valency + ')')
                name.append(anion[0][1])
                # formula = mp.reconstruct_formula(name.copy(), valency)
                hydr_i = material_string.find('hydrate')
                if hydr_i > -1:
                    pref = []
                    while material_string[hydr_i - 1] != ' ' and hydr_i > 0:
                        pref.append(material_string[hydr_i - 1])
                        hydr_i -= 1

                    pref = ''.join([p for p in reversed(pref)])

                    if pref not in self.__neg_prefixes:
                        # formula = formula+'·'+str(prefixes2num[pref])+'H2O'
                        name.append(pref + 'hydrate')
                result.append((''.join([n + ' ' for n in name]).strip(' '), valency))

        return result

    ###################################################################################################################
    # Misc
    ###################################################################################################################

    def cleanup_name(self, material_name):
        """
        cleaning up material name - fix due to tokenization imperfectness
        :param material_name: <str> material string
        :return: <str> updated material string
        """

        # print (material_name)

        material_name = material_name.replace('−', '-')

        if any(a in material_name for a in ['→', '⟶']):
            return ''

        if 'hbox' in material_name.lower():
            material_name = re.sub('(\\\\[a-z\(\)]+)', '', material_name)
            for t in ['{', '}', '_', ' ']:
                material_name = material_name.replace(t, '')
            material_name = material_name.rstrip('\\')

        material_name = re.sub('\({0,1}[0-9\.]*\s*⩽{0,1}\s*[x,y]{0,1}\s*[⩽=]\s*[0-9\.-]*\){0,1}', '', material_name)

        for c in ['\(99', '\(98', '\(90', '\(95', '\(96', '\(Alfa', '\(Aldrich', '\(A.R.', '\(Aladdin', '\(Sigma',
                  '\(A.G', '\(Fuchen', '\(Furuuchi', '\(AR\)']:
            split = re.split(c, material_name)
            if len(split) > 1:
                if len(split[1]) == '' or all(not s.isalpha() for s in split[1]):
                    material_name = re.split(c, material_name)[0]

        replace_dict = {'oxyde': 'oxide',
                        'luminum': 'luminium',
                        'magneshium': 'magnesium',
                        'stanate': 'stannate',
                        'sulph': 'sulf',
                        'buter': 'butyr',
                        'butir': 'butyr',
                        'butly': 'butyl',
                        'ethly': 'ethyl',
                        'ehtyl': 'ethyl',
                        'Abstract ': '',
                        'phio': 'thio',
                        'uim': 'ium',
                        'butryal': 'butyral',
                        'ooper': 'opper',
                        'acac': 'CH3COCHCOCH3',
                        'glass': '',
                        'glasses': '',
                        'europeam': 'europium',
                        'siliminite': 'sillimanite',
                        'acethylene': 'acetylene',
                        'iso-pro': 'isopro',
                        'anhydrous': ''
                        }

        for typo, correct in replace_dict.items():
            material_name = material_name.replace(typo, correct)

        if material_name[-2:] == '/C':
            material_name = material_name[:-2]

        material_name = re.sub('(poly[\s-])(?=[a-z])', 'poly', material_name)

        if material_name[-2:] == '/C':
            material_name = material_name[:-2]

        if len(material_name.split(' ')) > 1:
            for v in re.findall('[a-z](\([IV,]+\))', material_name):
                material_name = material_name.replace(v, ' ' + v)

        if material_name != '':
            material_name = self.__check_parentheses(material_name)

        material_name = material_name.lstrip(') -')
        material_name = material_name.rstrip('( ,.:;-±/δ')

        return material_name

    def __check_parentheses(self, formula):

        new_formula = formula

        new_formula = new_formula.replace('[', '(')
        new_formula = new_formula.replace(']', ')')
        new_formula = new_formula.replace('{', '(')
        new_formula = new_formula.replace('}', ')')

        par_open = re.findall('\(', new_formula)
        par_close = re.findall('\)', new_formula)

        if new_formula[0] == '(' and new_formula[-1] == ')' and len(par_close) == 1 and len(par_open) == 1:
            new_formula = new_formula[1:-1]

        if len(par_open) == 1 and len(par_close) == 0:
            if new_formula.find('(') == 0:
                new_formula = new_formula.replace('(', '')
            else:
                new_formula += ')'

        if len(par_close) == 1 and len(par_open) == 0:
            if new_formula[-1] == ')':
                new_formula = new_formula.rstrip(')')
            else:
                new_formula = '(' + new_formula
        #
        # if len(par_close) - len(par_open) == 1 and new_formula[-1] == ')':
        #     new_formula = new_formula.rstrip(')')

        return new_formula

    def __simplify(self, value):

        """
        simplifying stoichiometric expression
        :param value: string
        :return: string
        """

        for l in self.__greek_letters:
            _clash[l] = smp.Symbol(l)

        new_value = value
        for i, m in enumerate(re.finditer('(?<=[0-9])([a-z' + ''.join(self.__greek_letters) + '])', new_value)):
            new_value = new_value[0:m.start(1) + i] + '*' + new_value[m.start(1) + i:]
        new_value = smp.simplify(smp.sympify(new_value, _clash))
        if new_value.is_Number:
            new_value = round(float(new_value), 4)

        return str(new_value)

    def __lcm(self, x, y):
        """This function takes two
        integers and returns the L.C.M."""

        # choose the greater number
        lcm = None
        if x > y:
            greater = x
        else:
            greater = y

        found = False
        while not found:
            if (greater % x == 0) and (greater % y == 0):
                lcm = greater
                found = True
            greater += 1

        return lcm

    def __cast_stoichiometry(self, value):
        if value == 1:
            return ''

        return str(value)

    def __empty_structure(self):
        return dict(
            composition={},
            fraction_vars={},
            elements_vars={},
            hydrate='',
            phase=''
        )

    def __is_int(self, num):
        try:
            return round(float(num), 4) == round(float(num), 0)
        except:
            return False

    def __substitute_dopant(self, dopant, material_formula, material_composition):

        new_material_composition = {}
        new_material_formula = material_formula

        r = '[x0-9\.]+$'

        coeff = re.findall(r, dopant)
        element = re.split(r, dopant)[0]

        if coeff == [] or element not in self.__list_of_elements_1+self.__list_of_elements_2:
            return new_material_formula, material_composition

        for mat, compos in material_composition.items():
            expr = ''.join(['(' + v + ')+' for e, v in compos['composition'].items()]).rstrip('+')

            coeff = coeff[0] if not re.match('^[0]+[1-9]', coeff[0]) else '0.' + coeff[0][1:]
            expr = expr + '+(' + coeff + ')'
            #print ('-->', expr, self.__simplify(expr), self.__is_int(self.__simplify(expr)))

            if self.__is_int(self.__simplify(expr)):
                #print('-->', element, coeff)
                new_name = element + coeff + mat
                new_composition = compos['composition'].copy()
                new_composition[element] = coeff

                new_material_composition[new_name] = {}
                new_material_composition[new_name]['composition'] = new_composition
                new_material_composition[new_name]['fraction'] = compos['fraction']
                new_material_formula = new_material_formula.replace(mat, new_name)
            else:
                new_material_composition[mat] = {}
                new_material_composition[mat]['composition'] = compos['composition']
                new_material_composition[mat]['fraction'] = compos['fraction']

        return new_material_formula, new_material_composition