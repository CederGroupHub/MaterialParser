from text_cleanup import TextCleanUp
from materials_entity_recognition import MatRecognition
from material_parser.material_parser import MaterialParser
from chemdataextractor.doc import Paragraph
import re
import itertools

from pprint import pprint

class RecipeExtractor:
    def __init__(self, verbose=False, pubchem_lookup=False):
        self.__tp = TextCleanUp()
        self.__mer = MatRecognition()
        self.__mp = MaterialParser(pubchem_lookup=pubchem_lookup)
        self.__chemical_names = self.__mp.build_names_dictionary()
        self.__verbose = verbose
        self.__pubchem = pubchem_lookup

        self.__list_of_elements_1 = ['H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U']
        self.__list_of_elements_2 = ['He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'Cr',
                                     'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',
                                     'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'Xe',
                                     'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
                                     'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
                                     'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
                                     'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
                                     'Fl', 'Lv']

    def get_materials(self, doi, abstract_, syn_paragraph_, split_mixture=False):

        structures_list = []
        data_structure = {}

        # Text preprocessing: cleaning up paragraphs
        abstract = self.__tp.cleanup_text(abstract_)
        syn_paragraph = self.__tp.cleanup_text(syn_paragraph_)

        if self.__verbose:
            print('Abstract:', abstract)
            print('Synthesis paragraph:', syn_paragraph)

        #MER
        abstract_all, abstract_precursors, abstract_targets, abstract_others = self.__mer.mat_recognize(abstract)
        syn_all, syn_precursors, syn_targets, syn_others = self.__mer.mat_recognize(syn_paragraph)

        abstract_materials = list(set([m['text'] for m in abstract_all]))
        targets = list(set([m['text'] for m in syn_targets]))
        precursors = list(set([m['text'] for m in syn_precursors]))

        if self.__verbose:
            print('MER found materials')
            print('\tfrom abstract:', abstract_materials)
            print('\ttargets from synthesis paragaph', targets)
            print('\tprecursors from synthesis paragaph', precursors)

        #get dopants
        abstract_dopants = []
        abstract_materials_upd_1 = []
        for material in abstract_materials:
            dopants, new_name = self.__mp.get_dopants(material)

            abstract_dopants.extend(d for d in dopants)
            abstract_materials_upd_1.append(new_name)

        abstract_dopants = set(abstract_dopants)
        abstract_materials_upd_1 = list(set(abstract_materials_upd_1))

        targets_dopants = []
        targets_upd_1 = []
        for material in targets:
            dopants, new_name = self.__mp.get_dopants(material)

            targets_dopants.extend(d for d in dopants)
            targets_upd_1.append(new_name)

        targets_dopants = set(targets_dopants)
        targets_upd_1 = list(set(targets_upd_1))

        data_structure['dopants'] = list(targets_dopants.union(abstract_dopants))

        if self.__verbose:
            print('Removed doped parts:')
            print('\tdopants:', data_structure['dopants'])
            print('\tabstract materials', abstract_materials_upd_1)
            print('\ttargets', targets_upd_1)

        #build abbreviations dictionary
        materials = list(set(abstract_materials_upd_1+targets_upd_1))
        abbreviations = self.__mp.build_abbreviations_dict(materials, [abstract, syn_paragraph])

        if self.__verbose:
            print('Abbreviations dictionary:')
            pprint(abbreviations)

        abstract_materials_upd_2 = []
        for material in abstract_materials_upd_1:
            if material in abbreviations:
                abstract_materials_upd_2.append(abbreviations[material])
            else:
                abstract_materials_upd_2.append(material)

        targets_upd_2 = []
        for material in targets_upd_1:
            if material in abbreviations:
                targets_upd_2.append(abbreviations[material])
            else:
                targets_upd_2.append(material)

        abstract_materials_upd_2 = list(set(abstract_materials_upd_2))
        targets_upd_2 = list(set(targets_upd_2))

        if self.__verbose:
            print('After abbreviations substitution:')
            print('\tabstract:', abstract_materials_upd_2)
            print('\ttargets:', targets_upd_2)

        #checking valency
        valency = ''
        for material in abstract_materials_upd_2+targets_upd_2:
            v = re.findall('\s*\(([IV,]+)\)', material)
            if len(v) > 0:
                valency = v[0]

        #clean up materials names
        abstract_materials_upd_3 = []
        for material in abstract_materials_upd_2:
            new_name = self.__tp.clean_up_material_name(material)
            new_name = self.split_material_name(new_name)
            if new_name != '':
                abstract_materials_upd_3.append(new_name)
            else:
                abstract_materials_upd_3.append(material)
        abstract_materials_upd_3 = list(set(abstract_materials_upd_3))

        targets_upd_3 = []
        for material in targets_upd_2:
            new_name = self.__tp.clean_up_material_name(material)
            new_name = self.split_material_name(new_name)
            if new_name != '':
                targets_upd_3.append(new_name)
            else:
                targets_upd_3.append(material)
        targets_upd_3 = list(set(targets_upd_3))

        precursors_upd_3 = []
        for material in precursors:
            new_name = self.__tp.clean_up_material_name(material)
            new_name = self.split_material_name(new_name)
            if new_name != '':
                precursors_upd_3.append(new_name)
            else:
                precursors_upd_3.append(material)
        precursors_upd_3 = list(set(precursors_upd_3))

        if self.__verbose:
            print("After cleaning up names:")
            print ('\tabstract:', abstract_materials_upd_3)
            print ('\ttargets:', targets_upd_3)
            print ('\tprecursors:', precursors_upd_3)

        # parse chemical formulas of materials
        if self.__verbose:
            print('Extracting chemical composition:')

        abstract_materials_struct = {}
        for material in abstract_materials_upd_3:
            t_struct = {}
            try:
                t_struct = self.__mp.get_chemical_structure(material)
            except:
                t_struct['composition'] = {}
                print ('Failed: '+material)
            if len(t_struct['composition']) > 1:
                abstract_materials_struct[material] = t_struct

        targets_struct = {}
        for material in targets_upd_3:
            t_struct = {}
            try:
                t_struct = self.__mp.get_chemical_structure(material)
            except:
                t_struct['composition'] = {}
                print ('Failed: ' + material)
            if len(t_struct['composition']) > 1:
                targets_struct[material] = t_struct

        precursors_struct = {}
        for material in precursors_upd_3:
            t_struct = {}
            try:
                t_struct = self.__mp.get_chemical_structure(material)
            except:
                t_struct['composition'] = {}
                print ('Failed: ' + material)
            if len(t_struct['composition']) > 1 and \
                    all(el in self.__list_of_elements_1+self.__list_of_elements_2 for el in t_struct['composition']):
                precursors_struct[material] = t_struct

        #resolving fraction variables
        for material, struct in abstract_materials_struct.items():
            for var in struct['fraction_vars']:
                values, mode = self.get_values(var, syn_paragraph)
                if len(values) == 0:
                    values, mode = self.get_values(var, abstract)
                struct['fraction_vars'][var] = values

            for var in struct['elements_vars']:
                values, mode = self.get_values(var, syn_paragraph, elements=True)
                if len(values) == 0:
                    values, mode = self.get_values(var, abstract, elements=True)
                struct['elements_vars'][var] = values

        for material, struct in targets_struct.items():
            for var in struct['fraction_vars']:
                values, mode = self.get_values(var, syn_paragraph)
                if len(values) == 0:
                    values, mode = self.get_values(var, abstract)
                struct['fraction_vars'][var] = values

            for var in struct['elements_vars']:
                values, mode = self.get_values(var, syn_paragraph, elements=True)
                if len(values) == 0:
                    values, mode = self.get_values(var, abstract, elements=True)
                struct['elements_vars'][var] = values

        if self.__verbose:
            print('Materials structures from abstract:')
            pprint(abstract_materials_struct)
            print('Materials structures from targets:')
            pprint(targets_struct)
            print('Materials structures from precursors:')
            pprint(precursors_struct)

        #assigning possible abbreviations:
        for name, struct in abstract_materials_struct.items():
            struct['is_abbreviation_like'] = self.is_abbreviation_like(struct)

        for name, struct in targets_struct.items():
            struct['is_abbreviation_like'] = self.is_abbreviation_like(struct)

        final_targets = {name: struct for name, struct in targets_struct.items() if struct['composition'] != {} and not struct['is_abbreviation_like']}
        for name, struct in abstract_materials_struct.items():
            if not name in final_targets and struct['composition'] != {} and not struct['is_abbreviation_like']:
                final_targets[name] = struct

        #split mixtures if required
        if split_mixture:
            mixtures = [material for material, struct in final_targets.items() if struct['mixture'] != {}]
            for material in mixtures:
                struct = final_targets[material]
                for compound, compos in struct['mixture'].items():
                    if compound not in final_targets:
                        final_targets[compound] = {}
                        final_targets[compound]['chemical_name'] = ''
                        final_targets[compound]['composition'] = compos['composition']
                        final_targets[compound]['elements_vars'] = struct['elements_vars']
                        if any(not v.replace('.', '').replace('/','').isdigit() for v in compos['composition'].values()):
                            final_targets[compound]['fraction_vars'] = struct['fraction_vars']
                        else:
                            final_targets[compound]['fraction_vars'] = {}
                        final_targets[compound]['formula'] = compound
                        final_targets[compound]['name'] = struct['name']
                        final_targets[compound]['phase'] = struct['phase']
                        final_targets[compound]['mixture'] = {}
                del final_targets[material]


        if self.__verbose:
            print('List of final targets:')
            pprint(final_targets)

        for material, struct in final_targets.items():
            elements_array = self.get_substitutions_array(struct['elements_vars'])

            new_materials_array = []
            for subs in elements_array:
                composition = struct['composition'].copy()
                for var, val in subs.items():
                    composition[val] = composition[var]
                    del composition[var]
                new_materials_array.append(
                        dict(el_subst=subs, composition=composition, fraction_vars=struct['fraction_vars'].copy()))

            if len(new_materials_array) == 0:
                new_materials_array.append(dict(el_subst = {}, composition = struct['composition'].copy(), fraction_vars = struct['fraction_vars'].copy()))

            fractions_array = self.get_substitutions_array(struct['fraction_vars'])
            new_compositions = []
            for item in new_materials_array:
                if len(fractions_array) > 0:
                    for subs in fractions_array:
                        composition = item['composition'].copy()
                        obtained = True
                        for el, stoich in composition.items():
                            for var, val in subs.items():
                                stoich = stoich.replace(var, str(val))
                            try:
                                stoich = round(float(eval(stoich)), 3)
                                if stoich < 0:
                                    obtained = False
                            except:
                                obtained = False

                            composition[el] = str(stoich)

                        if obtained:
                            new_compositions.append(dict(
                                compos=composition,
                                subs=item['el_subst'],
                            ))

                else:
                    composition = item['composition'].copy()
                    obtained = True
                    for el, stoich in composition.items():
                        try:
                            stoich = round(float(eval(stoich)), 3)
                            if stoich < 0:
                                obtained = False
                        except:
                            obtained = False

                        composition[el] = str(stoich)

                    if obtained and all(composition != c['compos'] for c in new_compositions):
                        new_compositions.append(dict(
                            compos=composition,
                            subs=item['el_subst']
                        ))

            for composition in new_compositions:
                try:
                    correct_formula = self.reconstruct_formula(material, composition)
                except:
                    correct_formula = ''

                formula = ''
                zero_vals = []
                for el, val in composition['compos'].items():
                    if float(val) != 0.0:
                        formula = formula + el + self.__cast_stoichiometry(val)
                    else:
                        zero_vals.append(el)
                for el in zero_vals:
                    del composition['compos'][el]

                precursors = []
                for prec, p_struct in precursors_struct.items():
                    if self.__is_precursor(p_struct['composition'], composition['compos']):
                        precursors.append(prec)

                if len(precursors) > 1:
                    structures_list.append(dict(
                        init_name = material,
                        formula = formula,
                        correct_formula = correct_formula,
                        composition = composition['compos'],
                        substitutions = composition['subs'],
                        precursors = precursors,
                        precursors_compositions = {prec: precursors_struct[prec] for prec in precursors},
                        doi = doi,
                        syn_paragraph = syn_paragraph,
                        abstract = abstract
                    ))


        return structures_list

    def split_material_name(self, material_name, pubchem_lookup=False):

        """
        this is a fix of incorrect tokenization of chemical names
        do not rely much on it
        tested on specific sample of papers from 20K solid-state paragraphs
        use at your own risk
        :param
        material_name: string - initial material string
        pubchem_lookup: boolean - if True chemical name will be searched in pubchem DB
        :return: string
        """

        updated_name = ''
        # trying to split correctly name, irrelevant words and formula
        parts_1 = [w for w in re.split('\s', material_name) if len(w) > 2]
        parts_2 = [w for w in re.findall('[a-z]+', material_name) if len(w) > 2]

        if len(parts_1) > 1:
            if len(parts_1) == len(parts_2):
                if material_name.lower() in self.__chemical_names or material_name in self.__chemical_names:
                    updated_name = self.__chemical_names[material_name.lower()]
                # pubchem can be used to look for chemical names
                elif pubchem_lookup:
                    comp = pcp.get_compounds(material_name.lower(), 'name')
                    if len(comp) > 0:
                        updated_name = comp[0].molecular_formula

            if len(parts_1) - len(parts_2) == 1:
                try:
                    t_struct = self.get_structure_by_formula(parts_1[-1])
                    if t_struct['composition'] != {}:
                        updated_name = parts_1[-1]
                except:
                    updated_name = ''

            if len(updated_name) == 0:
                updated_name = material_name

        else:
            updated_name = material_name

        # in case there are parenthensis around
        if len(re.findall('[\[\]]', updated_name)) == 2 and updated_name[0] == '[' and updated_name[-1] == ']':
            updated_name = updated_name.replace('[', '')
            updated_name = updated_name.replace(']', '')

        return updated_name

    def get_values(self,var, text, elements=False):

        values = []
        mode = ''
        sents = Paragraph(text).sentences
        i = 0
        try:
            while len(values) == 0 and i < len(sents):
                if elements:
                    values = self.__mp.get_elements_values(var, sents[i].text.strip('., '))
                else:
                    values, mode = self.__mp.get_stoichiometric_values(var, sents[i].text.strip('.'))
                i = i + 1
        except:
            values = []
            mode = ''
            # fails_var.append(data)

        return values, mode

    def get_substitutions_array(self, subs_dict):

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

        if all(el[0].isupper() and len(el) == 1 and value == '1' for el, value in structure['composition'].items()) \
                and structure['composition'] != {}:
            return True

        if any(all(el[0].isupper() and len(el) == 1 and v == '1' for el, v in compos['composition'].items()) for compos in structure['mixture'].values()):
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

        if any(el in p_composition for el in t_composition.keys() if el != 'O') \
                and len(p_composition) <= len(t_composition) and p_composition != t_composition:
            return True

        return False

    def __cast_stoichiometry(self, value):

        value = float(value)
        if value == 1.0:
            return ''
        if value * 1000 % 1000 == 0.0:
            return str(int(value))

        return str(value)
