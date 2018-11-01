# coding=utf-8

__author__ = "Olga Kononova"
__maintainer__ = "Olga Kononova"
__email__ = "0lgaGkononova@yandex.ru"
__version__ = "1.0"

import regex as re
import collections
import sympy
from sympy.abc import _clash
import pubchempy as pcp
#from chemdataextractor.doc import Paragraph
import os
import json

# TODO
# Check if stoichiometry in () sums to 1.0 or integer

class MaterialParser:
    def __init__(self, pubchem_lookup=False):
        print('MaterialParser version 1.3')

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

        self.__chemical_names = self.build_names_dictionary()

        self.__ions = json.loads(open(os.path.join(self.__filename, 'rsc/ions_dictionary.json')).read())
        self.__anions = {ion['c_name']: {'valency': ion['valency'], 'e_name': ion['e_name'], 'n_atoms': ion['n_atoms']} for ion in self.__ions['anions']}
        self.__cations = {ion['c_name']: {'valency': ion['valency'], 'e_name': ion['e_name']} for ion in self.__ions['cations']}

        self.__prefixes2num = {'': 1, 'mono': 1, 'di': 2, 'tri': 3, 'tetra': 4, 'penta': 5, 'hexa': 6, 'hepta': 7,
                                'octa': 8,  'nona': 9, 'deca': 10 }

        self.__rome2num = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10}

        self.__pubchem_file = open(os.path.join(self.__filename, 'pubchem_log'), 'w')
        self.__pubchem_file.close()

        self.__pubchem = pubchem_lookup
        if pubchem_lookup:
            print ('Pubchem lookup is on! Will search for unknown materials name in PubChem DB.')

    ###################################################################################################################
    ### Methods to build chemical structure
    ###################################################################################################################

    def __simplify(self, value):

        """
        simplifying stoichiometric expression
        :param value: string
        :return: string
        """

        for l in self.__greek_letters:
            _clash[l] = sympy.Symbol(l)

        new_value = value
        for i, m in enumerate(re.finditer('(?<=[0-9])([a-z' + ''.join(self.__greek_letters) + '])', new_value)):
            new_value = new_value[0:m.start(1) + i] + '*' + new_value[m.start(1) + i:]
        new_value = sympy.simplify(sympy.sympify(new_value, _clash))
        if new_value.is_Number:
            new_value = round(float(new_value), 3)

        return str(new_value)

    def __get_sym_dict(self, f, factor):
        sym_dict = collections.defaultdict(str)
        r = "([A-Z]{1}[a-z]{0,1})\s*([-*\.\da-z" + ''.join(self.__greek_letters) + "\+/]*)"

        for m in re.finditer(r, f):
            """
            checking for correct elements names
            """
            el_bin = "{0}{1}".format(str(int(m.group(1)[0] in self.__list_of_elements_1 + ['M'])), str(
                int(m.group(1) in self.__list_of_elements_1 + self.__list_of_elements_2 + ['Ln', 'M'])))
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

    def __parse_parentheses(self, init_formula, init_factor, curr_dict):
        r = "\(((?>[^\(\)]+|(?R))*)\)\s*([-*\.\da-z\+/]*)"

        for m in regex.finditer(r, init_formula):
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
                new_formula = new_formula + ')'

        if len(par_close) == 1 and len(par_open) == 0:
            if new_formula[-1] == ')':
                new_formula = new_formula.rstrip(')')
            else:
                new_formula = '(' + new_formula

        if len(par_close) - len(par_open) == 1 and new_formula[-1] == ')':
            new_formula = new_formula.rstrip(')')

        return new_formula

    def get_structure_by_formula(self, formula):

        init_formula = formula
        formula = formula.replace(' ', '')
        formula = formula.replace('−', '-')

        formula = self.__check_parentheses(formula)

        elements_variables = collections.defaultdict(str)
        stoichiometry_variables = collections.defaultdict(str)

        # check for any weird syntax
        r = "\(([^\(\)]+)\)\s*([-*\.\da-z\+/]*)"
        for m in re.finditer(r, formula):
            if not m.group(1).isupper() and m.group(2) == '':
                formula = formula.replace('(' + m.group(1) + ')', m.group(1), 1)
            if ',' in m.group(1):
                elements_variables['M'] = re.split(',', m.group(1))
                formula = formula.replace('(' + m.group(1) + ')' + m.group(2), 'M' + m.group(2), 1)

        composition = self.__parse_formula(formula)

        # looking for variables in elements and stoichiometry
        for el, amt in composition.items():
            if el not in self.__list_of_elements_1 + self.__list_of_elements_2 + list(elements_variables.keys()):
                elements_variables[el] = []
            for var in re.findall('[a-z' + ''.join(self.__greek_letters) + ']', amt):
                stoichiometry_variables[var] = []

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

        # print("Check sum:")
        # to_calc = ''.join([composition[el] + '+' for el in composition])
        # to_calc = "(" + to_calc.rstrip('+- ') + ")"
        # print(self.__simplify(to_calc))

        formula_structure = dict(
            formula_=init_formula,
            formula=formula,
            composition=composition,
            fraction_vars=stoichiometry_variables,
            elements_vars=elements_variables,
        )

        return formula_structure

    def __empty_structure(self):
        return dict(
            composition={},
            mixture={},
            fraction_vars={},
            elements_vars={},
            formula='',
            chemical_name='',
            phase='',
            hydrate = ''
        )

    def __is_correct_composition(self, formula, chem_compos):
        if chem_compos == {}:
            return False
        if any(el not in formula + 'M' + 'Ln' + 'RE' + 'AE' or amt == '' for el, amt in chem_compos.items()):
            return False

        return True

    def __get_mixture(self, material_name):

        mixture = {}
        material_name = self.__check_parentheses(material_name)

        for m in re.finditer('\(1-[xy]{1}\)(.*)-\({0,1}[xy]{1}\){0,1}(.*)', material_name.replace(' ', '')):
            mixture[m.group(1)] = {}
            mixture[m.group(2)] = {}
            mixture[m.group(1)]['fraction'] = '1-x'
            mixture[m.group(1)]['composition'] = self.get_structure_by_formula(m.group(1))['composition']
            mixture[m.group(2)]['fraction'] = 'x'
            mixture[m.group(2)]['composition'] = self.get_structure_by_formula(m.group(2))['composition']

            for i in [1, 2]:
                if m.group(i)[0] == '(' and m.group(i)[-1] == ')':
                    line = m.group(i)[1:-1]
                    parts = [s for s in re.split('[-+]{1}([\d\.]*[A-Z][^-+]*)', line) if s != '' and s != line]
                    for s in parts:
                        name = re.findall('([\d\.]*)([\(A-Z].*)', s.strip(' -+'))[0]
                        mixture[name[1]] = {}
                        fraction = name[0]
                        if fraction == '': fraction = '1'
                        mixture[name[1]]['fraction'] = '(' + fraction + ')*' + mixture[m.group(i)]['fraction']
                        mixture[name[1]]['composition'] = self.get_structure_by_formula(name[1])['composition']
                    del mixture[m.group(i)]

        if mixture == {}:
            parts = [s for s in re.split('[-+]{1}([\d\.]*[A-Z][^-+]*)', material_name.replace(' ', '')) if
                     s != '' and s != material_name.replace(' ', '')]
            for s in parts:
                name = re.findall('([\d\.]*)([\(A-Z].*)', s.strip(' -+'))[0]
                mixture[name[1]] = {}
                fraction = name[0]
                if fraction == '': fraction = '1'
                mixture[name[1]]['fraction'] = fraction
                mixture[name[1]]['composition'] = self.get_structure_by_formula(name[1])['composition']

        for item in mixture:
            mixture[item]['fraction'] = self.__simplify(mixture[item]['fraction'])

        return mixture

    def get_chemical_structure(self, material_name):
        """
        The main function to obtain the closest chemical structure associated with a given material name
        :param material_name: string of material name
        :return: dictionary composition and stoichiometric variables
        """

        self.__pubchem_file = open(os.path.join(self.__filename, 'pubchem_log'), 'a')

        chemical_structure = self.__empty_structure().copy()

        # making sure the compound if not just ion
        if material_name[-1] == '+':
            return self.__empty_structure()

        # is there any phase specified
        phase = ''
        if material_name[0].islower():
            for m in re.finditer('([a-z' + ''.join(self.__greek_letters) + ']*)-{0,1}(.*)', material_name):
                if m.group(2) != '':
                    phase = m.group(1)
                    material_name = m.group(2)
                #print(m.group(1))
                #print(m.group(2))

        #print('Phase extraction: '+material_name)

        #material_name = re.sub('[' + chr(183) + '](.*)', '', material_name)

        # unifying hydrates names
        if 'hydrate' in material_name:
            for m in re.finditer('(.*) ([a-z]*)hydrate', material_name):
                if m.group(1) != '':
                    material_name = m.group(1)
                    chemical_structure['hydrate'] = self.__prefixes2num[m.group(2)]

        dots = [8901, 8729, 65381, 120, 42, 215, 8226]
        if 'H2O' in material_name:
            for c in dots:
                material_name = material_name.replace(chr(c), chr(183))
            for m in re.finditer('(.*)·([1-9]){0,1}H2O', material_name):
                if m.group(1) != '':
                    material_name = m.group(1)
                    chemical_structure['hydrate'] = m.group(2)

        # checking valency
        valency = [v.strip('()') for v in re.findall('(\s*\([IV,]+\))', material_name)]
        material_name = re.sub('(\s*\([IV,]+\))', '', material_name)

        #print('Symbol substitution: ' + material_name)

        # splitting mixture into compounds
        chemical_structure['mixture'] = self.__get_mixture(material_name)

        if chemical_structure['mixture'] != {}:
            chemical_structure['formula'] = ''.join(chemical_structure['mixture'].keys())
        else:
            chemical_structure['formula'] = material_name

        # trying to extract chemical structure for entire mixture or initial materials name
        try:
            t_struct = self.get_structure_by_formula(chemical_structure['formula'])
        except:
            t_struct = self.__empty_structure().copy()

        #if t_struct['composition'] != {}:
        chemical_structure['formula'] = t_struct['formula']
        chemical_structure['composition'] = t_struct['composition']
        chemical_structure['fraction_vars'] = t_struct['fraction_vars']
        chemical_structure['elements_vars'] = t_struct['elements_vars']

        # TODO: merge fraction variables

        # if material name is not proper formula look for it in dictionary
        # if not self.__is_correct_composition(chemical_structure['formula'], chemical_structure['composition']):
        #     if material_name.lower() in self.__chemical_names or material_name in self.__chemical_names:
        #         formula = self.__chemical_names[material_name.lower()]
        #         formula = re.sub('[' + chr(183) + '](.*)', '', formula)
        #         chemical_structure['elements_vars'] = collections.defaultdict(str)
        #         chemical_structure['fraction_vars'] = collections.defaultdict(str)
        #         try:
        #             t_struct = self.get_structure_by_formula(formula)
        #         except:
        #             t_struct = self.__empty_structure()
        #
        #         #if t_struct['composition'] != {}:
        #         chemical_structure['composition'] = t_struct['composition']
        #         chemical_structure['formula'] = t_struct['formula']
        #         chemical_structure['fraction_vars'] = t_struct['fraction_vars']
        #         chemical_structure['elements_vars'] = t_struct['elements_vars']

        # if material is not proper formula
        # first, trying to reconstruct formula and composition by name
        chem_terms = material_name.split(' ')
        if  len(chem_terms) == 1 and chem_terms[0].lower() in self.__cations|self.__anions:
            return self.__empty_structure()

        formula = self.reconstruct_by_name(chem_terms, valency)
        if formula != '':
            t_struct = self.get_structure_by_formula(formula)
            chemical_structure['formula'] = t_struct['formula']
            chemical_structure['composition'] = t_struct['composition']
            chemical_structure['fraction_vars'] = t_struct['fraction_vars']
            chemical_structure['elements_vars'] = t_struct['elements_vars']

        # second, trying PubChem
        if chemical_structure['composition'] == {} and self.__pubchem:
            chemical_structure['composition'] = collections.defaultdict(str)
            chemical_structure['elements_vars'] = collections.defaultdict(str)
            chemical_structure['fraction_vars'] = collections.defaultdict(str)
            chemical_structure['formula'] = ''

            pcp_compounds = pcp.get_compounds(material_name, 'name')
            if len(pcp_compounds) > 0:
                try:
                    formula = str(pcp_compounds[0].molecular_formula)
                    t_struct = self.get_structure_by_formula(formula)
                    # if len(t_struct['composition']) > 1:
                    #     print('Found in Pubchem:')
                    #     print(material_name)
                    #     print(pcp_compounds[0].molecular_formula)
                    self.__pubchem_file.write(str(material_name) + ' ' + formula+'\n')
                        #print (t_struct)
                except:
                    t_struct = self.__empty_structure().copy()

                chemical_structure['composition'] = t_struct['composition']
                chemical_structure['formula'] = t_struct['formula']

        chemical_structure['phase'] = phase

        # finally, check if there are variables in mixture
        for part in chemical_structure['mixture'].values():
            for var in re.findall('[a-z' + ''.join(self.__greek_letters) + ']', part['fraction']):
                chemical_structure['fraction_vars'][var] = []

        self.__pubchem_file.close()

        return chemical_structure

    # TODO method merging materials with same composition

    ###################################################################################################################
    ### Methods to substitute variables
    ###################################################################################################################

    def __get_values(self, string, mode, count=None, default_value=0.1, incr=None):
        values = []

        """
        given range
        """
        if mode == 'range' and len(string) != 0:
            string = string[0]
            if any(c in string for c in ['-', '–']):
                interval = re.split('[-–]', string)
            else:
                interval = [string[0].strip(' '), string[1].strip(' ')]

            if len(interval) > 0:
                values = [round(float(interval[0]), 4), round(float(interval[1]), 4)]
                # start = float(interval[0])
                # end = float(interval[1])
                # if incr != None:
                #     values = [round(start + i * incr, 4) for i in range(round((end - start) / incr))]
                #     values.append(interval[1])
                #
                # if count != None:
                #     incr = (end - start) / count
                #     values = [round(start + i * incr, 4) for i in range(count)]
                #     values.append(interval[1])
                #
                # if incr == None and count == None:
                #     values = [default_value]

        """
        given list
        """
        if mode == 'values' and len(string) != 0:
            values = [round(float(sympy.simplify(c)), 4) for c in re.findall('[0-9\./-]+', string[0].strip(' '))]

        return values

    def get_stoichiometric_values(self, var, sentence):
        values = []
        mode = ''
        # equal to exact values
        if len(values) == 0:
            values = self.__get_values(re.findall(var + '\s*=\s*([-]{0,1}[0-9\.\,/and\s]+)[\s\)\]\,]', sentence), mode='values')
            mode = 'values'
        # equal to range
        if len(values) == 0:
            values = self.__get_values(re.findall(var + '\s*=\s*([0-9-–\.\s]+)[\s\)\]\,m\%]', sentence), mode='range',
                                       count=5)
            mode = 'range'
        # within range
        if len(values) == 0:
            values = self.__get_values(
                re.findall('([0-9\.\s]*)\s*[<≤⩽]{0,1}\s*' + var + '\s*[<≤⩽>]{1}\s*([0-9.\s]+)[\s\)\]\.\,]', sentence),
                mode='range', count=5)
            mode = 'range'

        # from ... to ...
        if len(values) == 0:
            values = self.__get_values(
                re.findall(var + '[a-z\s]*from\s([0-9\./]+)\sto\s([0-9\./]+)', sentence), mode='range', count=5)
            mode = 'range'
        #
        # if values == []:
        #     values = [0.0]

        return values, mode

    def get_elements_values(self, var, sentence):
        values = re.findall(var + '\s*[=:]{1}\s*([A-Za-z,\s]+)', sentence)
        if len(values) > 0:
            values = [c for c in re.split('[,\s]', values[0]) if
                      c in self.__list_of_elements_1 + self.__list_of_elements_2]

        return values

    ###################################################################################################################
    ### Methods to resolve abbreviations
    ###################################################################################################################

    def __is_abbreviation(self, word):
        if all(c.isupper() for c in re.sub('[0-9x\-\(\)\.]', '', word)) and len(re.findall('[A-NP-Z]', word)) > 1:
            return True

        return False

    def build_abbreviations_dict(self, materials_list, paragraphs):
        """

        :param materials_list: list of found materials entities
        :param paragraphs: list of paragraphs where look for abbreviations
        :return: dictionary abbreviation - corresponding entity
        """

        abbreviations_dict = {t: '' for t in materials_list if self.__is_abbreviation(t.replace(' ', ''))}
        not_abbreviations = list(set(materials_list) - set(abbreviations_dict.keys()))

        # run through all materials list to resolve abbreviations among its entities
        for abbr in abbreviations_dict.keys():

            for material_name in not_abbreviations:
                if sorted(re.findall('[A-NP-Z]', abbr)) == sorted(re.findall('[A-NP-Z]', material_name)):
                    abbreviations_dict[abbr] = material_name

        # for all other abbreviations going through the paper text
        for abbr, name in abbreviations_dict.items():

            if name == '':

                sents = ' '.join(
                    [s.text for p in paragraphs for s in Paragraph(p).sentences if abbr in s.text]).split(abbr)
                i = 0
                while abbreviations_dict[abbr] == '' and i < len(sents):
                    sent = sents[i]
                    for tok in sent.split(' '):
                        if sorted(re.findall('[A-NP-Z]', tok)) == sorted(re.findall('[A-NP-Z]', abbr)):
                            abbreviations_dict[abbr] = tok
                    i = i + 1

        for abbr in abbreviations_dict.keys():
            parts = re.split('-', abbr)
            if all(p in abbreviations_dict for p in parts) and abbreviations_dict[abbr] == '' and len(parts) > 1:
                name = ''.join('('+abbreviations_dict[p]+')'+'-' for p in parts).rstrip('-')
                abbreviations_dict[abbr] = name


        empty_list = [abbr for abbr, name in abbreviations_dict.items() if name == '']
        for abbr in empty_list:
            del abbreviations_dict[abbr]

        return abbreviations_dict

    # def substitute_abbreviations(self, materials_list, abbreviations_dict):
    #
    #     updated_list = []
    #     for name in materials_list:
    #         if name in abbreviations_dict:
    #             updated_list.append(abbreviations_dict[name])
    #         else:
    #             updated_list.append(name)
    #
    #     return list(set(m for m in updated_list if m != ''))

    ###################################################################################################################
    ### Methods to separate doped part
    ###################################################################################################################

    def get_dopants(self, material_name):
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

        if any(w in new_material_name for w in ['mol', 'wt']):
            parts = re.split('[\-+]{0,1}[0-9\.]*[xyz]{0,1}\s*[mol]{3}\.{0,1}\%{0,1}', new_material_name) + \
                    re.split('[\-+]{0,1}[0-9\.]*[xyz]{0,1}\s*[wt]{2}\.{0,1}\%{0,1}', new_material_name)
            parts = sorted(list(set(parts) - set([new_material_name])), key=lambda x: len(x), reverse=True)

        if len(parts) > 1:
            new_material_name = parts[0].strip(' -+')
            dopants.extend(d.strip() for d in parts[1:] if d != '')

        # checking for dopants in formula
        parts = sorted(new_material_name.split(':'), key=lambda x: len(x), reverse=True)
        if len(parts) > 1:
            new_material_name = parts[0]
            dopants.append(parts[1])

        return dopants, new_material_name

    ###################################################################################################################
    ### Methods to call dictionary of inorganic compounds
    ###################################################################################################################

    def get_compounds_dictionary(self):
        names_dict = {}

        for line in open(os.path.join(self.__filename, 'rsc/inorganic_compounds_dictionary')):
            name, formula = line.strip().split(' – ')
            name = name[0].lower() + name[1:]
            names_dict[name] = formula

        for line in open(os.path.join(self.__filename, 'rsc/pub_chem_dictionary')):
            name, formula = line.strip().split(' – ')
            name = name[0].lower() + name[1:]
            names_dict[name] = formula

        formulas_dict = {formula: [] for formula in names_dict.values()}
        for name, formula in names_dict.items():
            formulas_dict[formula].append(name)

        return names_dict, formulas_dict

    def build_names_dictionary(self):

        names_dict, formulas_dict = self.get_compounds_dictionary()
        chemical_names = {}
        for name in names_dict.keys():
            chemical_names[name] = names_dict[name]
            chemical_names[re.sub('(\s*\([IV,]+\))', '', name)] = names_dict[name]

        return chemical_names

    def reconstruct_by_name(self, terms_list, name_valency_):

        chem_name = ''.join(c+' ' for c in terms_list).strip(' ')

        output_formula = ''

        if len(terms_list) < 2:
            return output_formula

        anion = terms_list.pop().lower()
        next_term = terms_list.pop().lower()
        guess_valency = 0
        name_valency = [self.__rome2num[v.upper()] for v in name_valency_]
        if len(name_valency) > 1:
            print ('WARNING! Found many valencies per chemical name '+ chem_name)
            print (name_valency)

        if 'hydrogen' in next_term:
            anion = next_term + ' ' + anion
        else:
            terms_list = terms_list + [next_term]
            for p in self.__prefixes2num.keys():
                if anion.find(p) == 0 and p != '':
                    guess_valency = self.__prefixes2num[p]
                    anion = re.sub('^'+p, '', anion).strip('-')
                    if anion == 'xide':
                        anion = 'o'+anion

        if anion in self.__anions:
           anion_data = self.__anions[anion]
        else:
            return output_formula

        if len(terms_list) != 1:
            return output_formula

        if terms_list[0] in self.__cations:
            cation_data = self.__cations[terms_list[0]]
        else:
            return output_formula

        true_valency = -1
        formulas = []
        for i, v_cation in enumerate(cation_data['valency']):
            v_anion = abs(anion_data['valency'][0])

            if v_cation in name_valency:
                true_valency = i

            cm = self.__lcm(v_cation, v_anion)
            v_cation = cm // v_cation
            v_anion = cm // v_anion

            if v_anion == guess_valency and guess_valency != 0:
                if true_valency not in [-1, i]:
                    print ('WARNING! Valencies mismatch for chemical name ' + chem_name)
                    print (name_valency, guess_valency)

                true_valency = i

            el_anion = anion_data['e_name']
            if v_anion > 1 and anion_data['n_atoms'] > 1:
                el_anion = '(' + el_anion + ')'

            formulas.append(cation_data['e_name'] + self.__cast_stoichiometry(v_cation)
                            + el_anion + self.__cast_stoichiometry(v_anion))

        true_valency = true_valency * int(len(cation_data['valency']) > 1)

        if true_valency == -1:
            print ('WARNING! Ambiguous valency for chemical name '+ chem_name)
            print ('Stoichiometry of compound will be undefined.')

            output_formula = cation_data['e_name'] + self.__cast_stoichiometry(abs(anion_data['valency'][0])) \
                             + el_anion + 'x'

        else:
            output_formula = formulas[true_valency]

        return output_formula

    def __valencies_combinations(self, v_list):
        output = []
        output_num = 1
        for l in v_list:
            output_num = output_num * len(l)
        output_size = len(v_list)
        for i, l in enumerate(v_list):
            period = output_num // len(l)
            output.append([el for i in range(period) for el in l])

        output = np.reshape(output, (output_num, output_size))
        return output

    def __lcm(self, x, y):
        """This function takes two
        integers and returns the L.C.M."""

        # choose the greater number
        if x > y:
            greater = x
        else:
            greater = y

        found = False
        while (not found):
            if ((greater % x == 0) and (greater % y == 0)):
                lcm = greater
                found = True
            greater += 1

        return lcm

    def __cast_stoichiometry(self, value):
        if value == 1:
            return ''

        return str(value)

    ###################################################################################################################
    ### Misc
    ###################################################################################################################







