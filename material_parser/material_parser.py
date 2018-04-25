import re
import regex
import collections
import sympy
from sympy.abc import _clash
import pubchempy as pcp
from synthesis_project_ceder.nlp.preprocessing import TextPreprocessor
from chemdataextractor.doc import Paragraph

import pprint

__author__ = "Olga Kononova"
__maintainer__ = "Olga Kononova"
__email__ = "0lgagkononova@gmail.com"


# TODO
# Rising errors
# Check if stoichiometry in () sums to 1.0 or integer
# Correct () pairs


class MaterialParser:
    def __init__(self):
        self.__list_of_elements_1 = ['H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U']
        self.__list_of_elements_2 = ['He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'Cr',
                                     'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',
                                     'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'Xe',
                                     'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
                                     'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
                                     'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
                                     'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
                                     'Fl', 'Lv']
        self.__list_of_trash_words = ['bulk', 'coated', 'rare', 'earth', 'ceramics', 'undoped']
        self.__greek_letters = ['α','δ', 'χ']

        self.__tp = TextPreprocessor('')

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
        for i, m in enumerate(re.finditer('(?<=[0-9])([a-z'+''.join(self.__greek_letters)+'])', new_value)):
            new_value = new_value[0:m.start(1) + i] + '*' + new_value[m.start(1) + i:]
        new_value = sympy.simplify(sympy.sympify(new_value, _clash))
        if new_value.is_Float:
            new_value = round(float(new_value), 3)

        return str(new_value)

    def get_sym_dict(self, f, factor):
        sym_dict = collections.defaultdict(str)
        r = "([A-Z]{1}[a-z]{0,1})\s*([-*\.\da-z"+''.join(self.__greek_letters)+"\+/]*)"

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

        unit_sym_dict = self.get_sym_dict(init_formula, init_factor)
        for el, amt in unit_sym_dict.items():
            if len(curr_dict[el]) == 0:
                curr_dict[el] = amt
            else:
                curr_dict[el] = '(' + str(curr_dict[el]) + ')' + '+' + '(' + str(amt) + ')'

        return curr_dict

    def parse_formula(self, init_formula):

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

        if new_formula[0] == '(' and new_formula[-1] == ')' and len(par_close) == 1 and  len(par_open) == 1:
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

        composition = self.parse_formula(formula)

        # looking for variables in elements and stoichiometry
        for el, amt in composition.items():
            if el not in self.__list_of_elements_1 + self.__list_of_elements_2 + list(elements_variables.keys()):
                elements_variables[el] = []
            for var in re.findall('[a-z'+''.join(self.__greek_letters)+']', amt):
                stoichiometry_variables[var] = []

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
            targets=[]
        )

        return formula_structure

    def __empty_structure(self):
        return dict(
            formula_='',
            formula='',
            composition=collections.defaultdict(str),
            stoichiometry_vars=collections.defaultdict(str),
            elements_vars=collections.defaultdict(str),
            targets=[]
        )

    def __is_correct_composition(self, formula, chem_compos):
        if chem_compos == {}:
            return False
        if any(el not in formula + 'M' or amt == '' for el, amt in chem_compos.items()):
            return False

        return True

    def get_mixture(self, material_name):

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
                        mixture[name[1]]['fraction'] = '('+fraction+')*'+mixture[m.group(i)]['fraction']
                        mixture[name[1]]['composition'] = self.get_structure_by_formula(name[1])['composition']
                    del mixture[m.group(i)]

        if mixture == {}:
            parts = [s for s in re.split('[-+]{1}([\d\.]*[A-Z][^-+]*)', material_name.replace(' ', '')) if s != '' and s != material_name.replace(' ', '')]
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
        chemical_structure = dict(
            composition={},
            mixture={},
            fraction_vars={},
            elements_vars={},
            formula='',
            chemical_name=''
        )

        material_name = re.sub('[∙⋅](.*)', '', material_name)

        chemical_structure['mixture'] = self.get_mixture(material_name)

        chemical_structure['formula'] = ''.join(chemical_structure['mixture'].keys())
        if chemical_structure['formula'] == '':
            chemical_structure['formula'] = material_name

        # trying to extract chemical structure
        try:
            t_struct = self.get_structure_by_formula(chemical_structure['formula'])
        except:
            t_struct = self.__empty_structure()
            #print('Something went wrong!' + material_name)
            # return self.__empty_structure()

        chemical_structure['composition'] = t_struct['composition']
        chemical_structure['fraction_vars'] = t_struct['fraction_vars']
        chemical_structure['elements_vars'] = t_struct['elements_vars']
        # TODO: merge fraction variables

        # if material name is not proper formula look for it in DB (pubchem, ICSD)
        if not self.__is_correct_composition(chemical_structure['formula'], chemical_structure['composition']):
            # print('Looking in pubchem ' + material_name)
            chemical_structure['composition'] = collections.defaultdict(str)
            # chemical_structure['stoichiometry_vars'] = collections.defaultdict(str)
            chemical_structure['elements_vars'] = collections.defaultdict(str)

            pcp_compounds = pcp.get_compounds(material_name, 'name')
            if len(pcp_compounds) > 0:
                try:
                    t_struct = self.get_structure_by_formula(pcp_compounds[0].molecular_formula)
                    chemical_structure['composition'] = t_struct['composition']
                    chemical_structure['fraction_vars'] = t_struct['fraction_vars']
                    chemical_structure['elements_vars'] = t_struct['elements_vars']
                except:
                    chemical_structure['composition'] = collections.defaultdict(str)

        # if still cannot find composition look word by word <- this is a quick fix due to wrong tokenization,
        # will be removed probably
        if chemical_structure['composition'] == {}:
            # print('Looking part by part ' + material_name)
            for word in material_name.split(' '):
                try:
                    t_struct = self.get_structure_by_formula(word)
                except:
                    t_struct = self.__empty_structure()
                    #print('Something went wrong!' + material_name)
                if self.__is_correct_composition(word, t_struct['composition']):
                    chemical_structure['composition'] = t_struct['composition']
                    chemical_structure['fraction_vars'] = t_struct['fraction_vars']
                    chemical_structure['elements_vars'] = t_struct['elements_vars']

        chemical_structure['name'] = material_name

        # finally, check if there are variables in mixture
        for part  in chemical_structure['mixture'].values():
            for var in re.findall('[a-z'+''.join(self.__greek_letters)+']', part['fraction']):
                chemical_structure['fraction_vars'][var] = []

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
                start = float(interval[0])
                end = float(interval[1])
                if incr != None:
                    values = [round(start + i * incr, 4) for i in range(round((end - start) / incr))]
                    values.append(interval[1])

                if count != None:
                    incr = (end - start) / count
                    values = [round(start + i * incr, 4) for i in range(count)]
                    values.append(interval[1])

                if incr == None and count == None:
                    values = [default_value]

        """
        given list
        """
        if mode == 'values' and len(string) != 0:
            values = [round(float(sympy.simplify(c)), 4) for c in re.findall('[0-9\./]+', string[0].strip(' '))]

        return values

    def get_stoichiometric_values(self, var, sentence):
        values = []
        # equal to exact values
        if len(values) == 0:
            values = self.__get_values(re.findall(var + '\s*=\s*([0-9\.\,/and\s]+)[\s\)\]\,]', sentence), mode='values')
        # equal to range
        if len(values) == 0:
            values = self.__get_values(re.findall(var + '\s*=\s*([0-9-–\.\s]+)[\s\)\]\,m\%]', sentence), mode='range',
                                       count=5)
        # within range
        if len(values) == 0:
            values = self.__get_values(
                re.findall('([0-9\.\s]*)\s*[<≤⩽]{0,1}\s*' + var + '\s*[<≤⩽>]{1}\s*([0-9.\s]+)[\s\)\]\.\,]', sentence),
                mode='range', count=5)

        # from ... to ...
        if len(values) == 0:
            values = self.__get_values(
                re.findall(var + '[a-z\s]*from\s([0-9\./]+)\sto\s([0-9\./]+)', sentence), mode='range', count=5)

        return values

    def get_elements_values(self, var, sentence):
        values = re.findall(var + '\s*[=:]{1}\s*([A-Za-z,\s]+)', sentence)
        if len(values) > 0:
            values = [c for c in re.split('[,\s]', values[0]) if
                      c in self.__list_of_elements_1 + self.__list_of_elements_2]

        return values

    # def replace_stoichiometry(self, string, variable, amount):
    #     new_string = '0'
    #     if len(re.findall('([0-9]*[a-z]*)' + variable, string)) != 0:
    #         if re.findall('([0-9]*[a-z]*)' + variable, string)[0] != '':
    #             new_string = string.replace(variable, '*' + str(amount))
    #     else:
    #         new_string = string.replace(variable, str(amount))
    #
    #     return new_string
    #
    # def substitute_elements(self, structure):
    #     output = []
    #     for elem, values in structure['elements_vars'].items():
    #         # print('For element', elem, values)
    #         for v in values:
    #             target = dict(structure['composition'])
    #             target[v] = target[elem]
    #             del target[elem]
    #             output.append(dict(
    #                 formula_=structure['formula_'],
    #                 formula=structure['formula'],
    #                 composition=dict(target),
    #                 stoichiometry_vars=structure['stoichiometry_vars'],
    #                 elements_vars=structure['elements_vars'],
    #                 targets=structure['targets']
    #             ))
    #
    #     return output
    #
    # def substitute_stoichiometry(self, structure):
    #     output = []
    #     for var, amounts in structure['stoichiometry_vars'].items():
    #         for amt in amounts:
    #             target_var = dict(structure['composition'])
    #             for el in target_var:
    #                 target_var[el] = eval(self.replace_stoichiometry(target_var[el], var, amt))
    #             output.append(target_var)
    #
    #     return output
    #
    # def test_parsing(self, material_name, sentence):
    #     formula_structure = self.get_structure_by_formula(material_name)
    #
    #     # find stoichiometric variables in text
    #     for var in list(formula_structure['stoichiometry_vars'].keys()):
    #         if len(formula_structure['stoichiometry_vars'][var]) == 0:
    #             formula_structure['stoichiometry_vars'][var] = self.get_stoichiometric_values(var, sentence)
    #
    #     # find element variables
    #     for var in list(formula_structure['elements_vars'].keys()):
    #         if len(formula_structure['elements_vars'][var]) == 0:
    #             formula_structure['elements_vars'][var] = self.get_elements_values(var, sentence)
    #
    #     # substitute
    #     targets = self.substitute_elements(formula_structure)
    #
    #     if len(targets) == 0:
    #         targets.append(dict(
    #             formula_=formula_structure['formula_'],
    #             formula=formula_structure['formula'],
    #             composition=formula_structure['composition'],
    #             stoichiometry_vars=formula_structure['stoichiometry_vars'],
    #             elements_vars=formula_structure['elements_vars'],
    #             targets=formula_structure['targets']
    #         ))
    #
    #     for target in targets:
    #         formula_structure['targets'].extend(dict(t) for t in self.substitute_stoichiometry(target))
    #
    #     return formula_structure

    ###################################################################################################################
    ### Methods to resolve abbreviations
    ###################################################################################################################

    def __is_abbreviation(self, word):
        if all(c.isupper() for c in re.sub('[0-9x\-\(\)\.]', '', word)) and len(re.findall('[A-NP-Z]', word)) > 1:
            return True

        return False

    def get_abbreviations_dict(self, materials_list, paragraphs):
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

            abbr_ = self.__tp._preprocess_text(abbr)

            if name == '':

                sents = ' '.join(
                    [s.text for p in paragraphs for s in Paragraph(p).sentences if abbr_ in s.text]).split(abbr_)
                i = 0
                while abbreviations_dict[abbr] == '' and i < len(sents):
                    sent = sents[i]
                    for tok in sent.split(' '):
                        if sorted(re.findall('[A-NP-Z]', tok)) == sorted(re.findall('[A-NP-Z]', abbr_)):
                            abbreviations_dict[abbr] = tok
                    i = i + 1

        return abbreviations_dict

    def substitute_abbreviations(self, materials_list, abbreviations_dict):

        updated_list = []
        for name in materials_list:
            if name in abbreviations_dict:
                updated_list.append(abbreviations_dict[name])
            else:
                updated_list.append(name)

        return list(set(m for m in updated_list if m != ''))

    ###################################################################################################################
    ### Methods to separate doped part
    ###################################################################################################################

    def get_dopants(self, material_name):
        new_material_name = material_name
        dopants = []

        s = re.split('doped with', new_material_name)
        new_material_name = s[0].strip(' -+')
        if len(s) > 1:
            dopants.append(s[1].strip())

        # checking for doped element prefix
        new_material_name = new_material_name.replace('codoped', 'doped')
        dopant = ''
        r = "(.*)\s*[-]*\s*doped\s*"
        # r = "(\w*)-doped\s*"
        for m in re.finditer(r, new_material_name):
            dopant = m.group(1).strip(' -')
            new_material_name = new_material_name.replace(m.group(0), '')

        if dopant != '':
            dopants.append(dopant)

        s = re.split('[\-+]{0,1}[0-9\.]*[xyz]{0,1}\s*[mol]{3}\.{0,1}\%{0,1}', new_material_name) + \
            re.split('[\-+]{0,1}[0-9\.]*[xyz]{0,1}\s*[wt]{2}\.{0,1}\%{0,1}', new_material_name)
        s = sorted(list(set(s) - set([new_material_name])), key=lambda x: len(x), reverse=True)

        if len(s) > 0:
            new_material_name = s[0].strip(' -+')
            dopants.extend(d.strip() for d in s[1:] if d != '')

        # checking for dopants in formula
        parts = sorted(new_material_name.split(':'), key=lambda x: len(x), reverse=True)
        new_material_name = parts[0]
        if len(parts) > 1:
            dopants.append(parts[1])

        return dopants, new_material_name
