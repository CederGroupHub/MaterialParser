# coding=utf-8
import regex as re
import material_parser.core.regex_parser as rp
import material_parser.core.chemical_sets as cs
from material_parser.core.preprocessing_tools.preprocessing_abc import PreprocessingABC
from material_parser.core.formula_processing import process_formula
from material_parser.core.utils import lcm, cast_stoichiometry


class ChemicalNameProcessing(PreprocessingABC):

    def __init__(self, regex_parser):
        super(ChemicalNameProcessing, self).__init__(regex_parser)

    def process_string(self, material_string, chemical_structure):
        """
        attempts to reconstruct chemical formula from arbitrary chemical name
        :param material_string:
        :param chemical_structure:
        :return:
        """
        if not material_string:
            return material_string, chemical_structure

        """
        checking chemical name is default dictionary
        """
        material_string_nodashes = " ".join([t for t in re.split(r'\s|(?<=[a-z])-', material_string)])  # ???
        material_string_nodashes = material_string_nodashes[0].lower() + material_string_nodashes[1:]
        if material_string_nodashes in cs.pubchem_dictionary:
            chemical_structure.material_name = material_string_nodashes
            chemical_structure.material_formula = cs.pubchem_dictionary[material_string_nodashes]
            return material_string, chemical_structure

        """
        checking abbreviations and acronyms dictionary
        """
        if material_string in cs.default_abbreviations:
            chemical_structure.material_formula = cs.default_abbreviations[material_string]
            return material_string, chemical_structure

        """
        separate formula from the string if possible, e.g. praseodymium tungstate Cd0.25Pr0.50□0.25WO4
        """
        material_formula, material_name = self.__separate_formula_from_string(material_string)

        """
        attempt to reconstruct material formula from the string of terms
        """
        if material_formula == "":
            material_formula = self.__reconstruct_formula_from_string(material_name)

        material_name = "" if material_formula == "" else material_name
        material_string = material_formula if material_formula else material_string
        #material_formula = material_string if material_formula == "" else material_formula

        chemical_structure.material_name = material_name
        chemical_structure.material_formula = material_formula

        return material_string, chemical_structure

    def __separate_formula_from_string(self, material_string):
        """
        Separate chemical formula if presented among chemical terms.
        This happens due to poor tokenization
        :param material_string
        :return: formula: <str> chemical formula found in material string
                 name: <str> chemical name found in material string
        """
        material_string = material_string.replace("[", " [").strip()
        material_string = re.sub(r"\s{2,}", " ", material_string)
        split = re.split(r"\s", material_string)
        """
        if it is one term - return the term only
        """
        if len(split) == 1:
            return "", material_string

        """
        if two terms - checking if they are proper cation + anion
        """
        if len(split) == 2 \
            and all(t in cs.list_of_elements or t.rstrip('s').lower() in cs.list_of_anions|{'metal'} for t in split):
            t1 = cs.element2name[split[0]] if split[0] in cs.element2name else split[0].rstrip('s').lower()
            t2 = cs.element2name[split[1]] if split[1] in cs.element2name else split[1].rstrip('s').lower()
            return "", t1 + " " + t2

        """
        trying each part if it is a proper formula term
        """
        formulas = []
        terms = []
        for part in split:
            composition_part = self.__try_formula(part)
            if composition_part != {} or part.strip('+-1234567890') in cs.list_of_elements:
                formulas.append(part)
            else:
                terms.append(part)

        """
        more than one formula is found
        """
        if len(formulas) > 1:
            return material_string, ""

        """
        no any formula is found
        """
        if len(formulas) == 0:
            return "", material_string

        """
        string can be a formula split into pieces due to tokenization errors
        """
        if any(t.strip('+-1234567890') in cs.list_of_elements for t in terms):
            return material_string, ""

        """
        correctly found one formula and proper chemical terms
        """
        return formulas[0], " ".join(terms)

    def __reconstruct_formula_from_string(self, material_name, valency=""):
        """
        reconstructing chemical formula for simple chemical names anion + cation
        :param material_name: <str> chemical name
        :param valency: <str> anion valency
        :return: <str> chemical formula
        """

        output_formula = ""
        material_name = re.sub(r"(\([IV]*\))", " \\1 ", material_name)  # removing valence symbols
        material_name = re.sub(r"\s{2,}", " ", material_name)

        terms_list, valency_list = [], []
        hydrate = ""
        cation_prefix_num = 0
        cation_data = {"c_name": "",
                       "valency": [],
                       "e_name": "",
                       "n_atoms": 0}

        """
        checking if any of the terms represents valency or refers to hydrate
        """
        for t in material_name.split(" "):
            if t.strip("()") in cs.rome2num:
                valency_list.append(cs.rome2num[t.strip("()")])
                continue
            if "hydrate" in t.lower():
                hydrate = t
                continue
            terms_list.append(t.strip(" -"))

        """
        discarding prefixes. e.g. sulfur di-oxide
        """
        terms_list_upd = []
        for t in terms_list:
            if all(p not in cs.prefixes2num for p in t.split("-")):
                terms_list_upd.extend(_ for _ in t.split("-"))
            else:
                terms_list_upd.append("".join(t.split("-")))
        terms_list = terms_list_upd

        """
        if chemical term is the name of element/ion, return the element/ion
        """
        t = " ".join(terms_list).lower()
        if t in cs.anions:
            return cs.anions[t]["e_name"]
        if t in cs.cations:
            return cs.cations[t]["e_name"]

        """
        if only one term, no formula is available
        """
        if len(terms_list) < 2:
            return output_formula

        """
        otherwise, getting anion term
        """
        anion = terms_list.pop().lower().rstrip("s")

        """
        guessing valency
        """
        if valency == "":
            valency_num = max(valency_list + [0])
        else:
            valency_num = cs.rome2num[valency.strip("()")]

        next_term = terms_list.pop()
        """
        checking is the next term is referring to hydrogen, e.g. hydrogen phosphate
        """
        if "hydrogen" in next_term.lower() and len(terms_list) != 0:
            anion = next_term + " " + anion
        else:
            terms_list += [next_term]

        """
        checking for prefix in anion, e.g. dihydrogen phosphate
        """
        anion_prefix_num, anion = (0, anion) if anion.lower() in cs.anions else self.__get_prefix(anion)

        """
        getting chemical data for anion from ion_dictionary.json
        """
        if anion in cs.anions:
            anion_data = cs.anions[anion].copy()
        elif anion in ["metal"] and terms_list[0].lower() in cs.cations:
            return cs.cations[terms_list[0].lower()]["e_name"]
        else:
            return output_formula

        """
        if more than one term left - return empty string
        """
        if len(terms_list) >= 2:
            return output_formula

        """
        getting chemical data for cation from ion_dictionary.json
        """
        if len(terms_list) == 1:
            cation = terms_list[0]
            cation_prefix_num, cation = (0, cation) if cation.lower() in cs.cations else self.__get_prefix(cation)
            if cation.lower() in cs.cations:
                cation = cation.lower()
                cation_data = cs.cations[cation].copy()
            elif cation in cs.element2name:
                cation_data = cs.cations[cs.element2name[cation]].copy()
            else:
                return output_formula

        """
        print warning if found valency is not in list of available cation valencies
        """
        if len(cation_data["valency"]) > 1 and valency_num != 0:
            if valency_num not in cation_data["valency"]:
                print("WARNING! Not common valency value for " + material_name)
                print(cation_data["valency"])
                print(valency_num)
            cation_data["valency"] = [valency_num]

        """
        reconstructing formula
        """
        output_formula = self.__build_formula(anion=anion_data,
                                              cation=cation_data,
                                              cation_prefix_num=cation_prefix_num,
                                              anion_prefix_num=anion_prefix_num)

        """
        adding water molecule to hydrates
        """
        if hydrate != "":
            hydrate_prefix_num, hydrate = self.__get_prefix(hydrate)
            hydrate_prefix = "" if hydrate_prefix_num in [0, 1] else str(hydrate_prefix_num)
            output_formula = output_formula + "·" + hydrate_prefix + "H2O"

        return output_formula

    def __try_formula(self, formula):
        if re.match(r"(\s*\([IV,]+\))", formula):
            return {}

        #try:
        chem_struct = process_formula(formula, self._re)
        composition = {"elements": chem_struct["elements"],
                       "elements_x": chem_struct["elements_x"]}
        # except:
        #     composition = {'elements': {}, 'elements_x': {}}
        if self.__is_abbreviation_composition(composition):
            return {}

        return composition["elements"]

    @staticmethod
    def __is_abbreviation_composition(composition):
        if all(e.isupper() and s in ["1.0", "1", "s"] for e, s in composition["elements"].items()):
            return True

        elements_vars = [el for el in composition["elements_x"].keys() if len(el) == 1 and el.isupper()]
        if len(elements_vars) > 1:
            return True

        return False

    @staticmethod
    def __get_prefix(material_name):

        pref_num = 0
        material_name_upd = material_name

        for p in cs.prefixes2num.keys():
            if material_name.lower().find(p) == 0 and p != "":
                pref_num = cs.prefixes2num[p]
                material_name_upd = material_name_upd[len(p):].strip("-")

                if material_name_upd == "xide":
                    material_name_upd = "oxide"

        return pref_num, material_name_upd

    @staticmethod
    def __build_formula(cation, anion, cation_prefix_num=0, anion_prefix_num=0):
        """
        naive implementation of formula reconstruction from cation and anion data
        :param cation: {"c_name" : "sodium", "valency" : [1], "e_name" : "Na", "n_atoms" : 1}
        :param anion: {"c_name" : "oxide", "valency" : [-2], "e_name" : "O", "n_atoms" : 1}
        :param cation_prefix_num:
        :param anion_prefix_num:
        :return:
        """

        cation_stoich, anion_stoich = 0, 0

        if anion_prefix_num + cation_prefix_num == 0 or anion_prefix_num * cation_prefix_num != 0:
            v_cation, v_anion = abs(cation["valency"][0]), abs(anion["valency"][0])
            cm = lcm(v_cation, v_anion)
            cation_stoich, anion_stoich = cm // v_cation, cm // v_anion

        if anion_prefix_num != 0:
            cation_stoich, anion_stoich = 0, anion_prefix_num
            i = 0
            while cation_stoich == 0 and i < len(cation["valency"]):
                cation_stoich = anion_prefix_num * abs(anion["valency"][0]) // abs(cation["valency"][i])
                i = i + 1

        if cation_prefix_num != 0:
            cation_stoich = cation_prefix_num
            anion_stoich = cation_prefix_num * abs(cation["valency"][0]) // abs(anion["valency"][0])

        anion_name_el = anion["e_name"]
        if anion_stoich > 1 and anion["n_atoms"] > 1:
            anion_name_el = "(" + anion_name_el + ")"

        cation_name_el = cation["e_name"]
        if cation_stoich > 1 and cation["n_atoms"] > 1:
            cation_name_el = "(" + cation_name_el + ")"

        return "{0}{1}{2}{3}".format(cation_name_el,
                                     cast_stoichiometry(cation_stoich),
                                     anion_name_el,
                                     cast_stoichiometry(anion_stoich))
