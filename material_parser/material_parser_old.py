# coding=utf-8

__author__ = "Olga Kononova"
__maintainer__ = "Olga Kononova"
__email__ = "0lgaGkononova@yandex.ru"
__version__ = "7.0.0"

import material_parser.core.chemical_sets as cs
import material_parser.core.constants as cnst
import collections
import itertools
import json
import os
import pubchempy as pcp
import regex as re
import sympy as smp
from pprint import pprint
from sympy.abc import _clash

from material_parser.core.string_cleanup import cleanup_name, check_parentheses


class MaterialParser:
    def __init__(self, verbose=False, pubchem_lookup=False, fails_log=False, dictionary_update=False):

        self.__filename = os.path.dirname(os.path.realpath(__file__))

        # self.__fails_log = fails_log
        # if fails_log:
        #     self.__pubchem_file = open(os.path.join(self.__filename, "pubchem_log"), "w")
        #     self.__pubchem_file.close()

        self.__diatomic_molecules = {"O2": collections.OrderedDict([("O", "2")]),
                                   "N2": collections.OrderedDict([("N", "2")]),
                                   "H2": collections.OrderedDict([("H", "2")]),
                                   "H2O": collections.OrderedDict([("H", "2"), ("O", "1")])}

        self.__pubchem = pubchem_lookup
        if pubchem_lookup:
            print("Pubchem lookup is on! Will search for unknown materials name in PubChem DB.")

        self.__dictionary_update = dictionary_update
        if dictionary_update:
            self.__dictionary_file = open(os.path.join(self.__filename, "dictionary_update"), "w")
            self.__dictionary_file.close()

        self.__verbose = verbose

    ###################################################################################################################
    # Parsing material name
    ###################################################################################################################

    def parse_material_string(self, material_string_):
        """
        Main method to convert chemical name into chemical formula and parse formula into chemical structure
        :param material_string_: < str> - material name and/or formula
        :return: dict(material_string: <str> - initial material string,
                      material_name: <str> - chemical name of material if found in the string
                      material_formula: <str> - chemical formula of material reconstructed from the string
                      composition: <list> of Object(dict):
                                formula: <str> - compound formula
                                amount: <str> - fraction of compound in mixture
                                elements: Object(dict) - {element: ratio/amount}
                      additives: <list> of <str> - list of dopped materials/elements found in material string
                      amount_vars: Object(dict) - {amount variable: <list> - values}
                      elements_vars: Object(dict) - {element variable: <list> - values}
                      oxygen_deficiency: <str> - sign of fraction (excess or deficiency) as given in material string
                      phase: <str> - material phase given in material string
                      is_acronym: <bool> - material string looks like acronym
        """

        material_string = cleanup_name(material_string_)
        if self.__verbose:
            print("After cleaning up string:", material_string_, "-->", material_string)

        if material_string in cs.list_of_elements:
            return self.__element_structure(material_string)

        if material_string in self.__diatomic_molecules:
            output_structure = self.__element_structure(material_string)
            output_structure["composition"][0]["elements"] = self.__diatomic_molecules.get(material_string)
            return output_structure

        """
        Converting material string into chemical formula
        """
        if material_string in cs.name2element:
            return self.__element_structure(cs.name2element.get(material_string))

        phase, material_string = self.separate_phase(material_string)
        additives, material_string = self.separate_additives(material_string)
        if self.__verbose:
            print("After additives extraction:", material_string, "|", additives)

        _, material_formula, material_name = self.string2formula(material_string)

        """
        Extracting composition from chemical formula
        """
        material_compounds = self.split_formula_into_compounds(material_formula)

        if self.__verbose:
            print("\tAfter splitting:", material_formula, "-->", material_compounds)

        output_structure = self.__empty_structure().copy()
        output_structure['material_string'] = material_string_
        output_structure['material_name'] = material_name
        output_structure["material_formula"] = material_formula
        hydrate = [f for m, f in material_compounds if m == "H2O"]
        material_compounds = [(m, f) for m, f in material_compounds if m != "H2O"]
        oxygen_deficiency = None
        for compound, amount in material_compounds:
            if compound in cs.default_abbreviations:
                if self.__verbose:
                    print("Found abbreviation:", compound, "-->", cs.default_abbreviations.get(compound))
                compound = cs.default_abbreviations.get(compound)
            composition = self.formula2composition(compound)
            if self.__verbose:
                print("Resolved composition:", compound, "-->", composition)
            output_structure["amounts_vars"].update(composition.get("amounts_vars"))
            output_structure["elements_vars"].update(composition.get("elements_vars"))
            if composition["elements"] != {}:
                output_structure["composition"].append({"formula": composition.get("formula"),
                                                        "amount": amount,
                                                        "elements": composition.get("elements")})
            if composition["oxygen_deficiency"]:
                oxygen_deficiency = composition["oxygen_deficiency"]
            output_structure["phase"] = phase
        output_structure["oxygen_deficiency"] = oxygen_deficiency

        """
        assigning abbreviations
        """
        output_structure["is_acronym"] = self.__is_acronym(output_structure)

        """
        substituting additive into composition if it makes fractions to sum-up to integer
        adding to composition if it is compound
        """
        if len(additives) == 1 and not output_structure["is_acronym"]:
            output_structure = self.substitute_additives(additives, output_structure)
        else:
            output_structure["additives"] = additives

        """
        appending hydrate part
        """
        if hydrate != []:
            output_structure["composition"].append(
                {"formula": "H2O",
                 "amount": hydrate[0],
                 "elements": self.__diatomic_molecules["H2O"],
                 "species": self.get_species("H2O")
                 }
            )

        """
        correction for negative stoichiometry and amount
        """
        try:
            if any(float(s) < 0.0 or float(s) > 100.0
                   for compound in output_structure["composition"] for e, s in compound["elements"].items()):
                output_structure["composition"] = []
        except:
            pass
        for compound in output_structure["composition"]:
            if len(re.findall("[b-mo-w]+", compound["amount"])) > 0:
                compound["amount"] = "1"

        """
        extracting species
        """
        if not output_structure["is_acronym"]:
            for compound in output_structure["composition"]:
                try:
                    compound["species"] = self.get_species(compound["formula"]) \
                        if len(compound["elements"]) > 2 or compound["formula"] == "H2O" else compound["elements"]
                except:
                    compound["species"] = collections.OrderedDict()
        else:
            for compound in output_structure["composition"]:
                compound["elements"] = {}
                compound["species"] = collections.OrderedDict()
            output_structure["elements_vars"] = {}

        """
        finally, combine unified material formula from its compounds 
        """
        output_structure["material_formula"] = self.__combine_formula(output_structure["composition"])

        # if output_structure["composition"] == [] and self.__fails_log:
        #     with open(os.path.join(self.__filename, "fails_log"), "a") as f_log:
        #         f_log.write(material_name + "\n")

        return output_structure

    def formula2composition(self, formula):
        """
        Parsing chemical formula in composition
        :param formula: <str> chemical formula
        :return: dict(formula: <str> formula string corresponding to obtained composition
                     composition: <dict> element: fraction
                     fraction_vars: <dict> elements fraction variables: <list> values
                     elements_vars: <dict> elements variables: <list> values
                     hydrate: <str> if material is hydrate fraction of H2O
                     phase: <str> material phase appeared in formula
                    )
        """

        formula = formula.replace(" ", "")
        """
        find oxygen deficiency
        """
        formula, oxygen_deficiency, oxygen_deficiency_sym = self.__get_oxygen_deficiency(formula)
        """
        converting fractions a(b+x)/c into (a/c*b+a/c*x)
        """
        formula = self.__fraction_convertion(formula)

        elements_variables = collections.defaultdict(str)
        stoichiometry_variables = collections.defaultdict(str)
        """
        check for any weird syntax (A,B)zElxEly...
        replacing with MzElxEly... and M = [A, B]
        """
        r = r"(\([A-Za-z\s]+[\/,\s]+[A-Za-z]+\))"
        for m in re.finditer(r, formula):
            elements_variables["M"] = re.split(r"[\/,]", m.group(0).strip('()'))
            formula = formula.replace(m.group(0), "M", 1)

        composition = self.__parse_formula(formula)

        if re.findall("[a-z]{4,}", formula) != [] and composition != {}:
            composition = collections.OrderedDict()
        """
        looking for variables in elements and stoichiometry 
        """
        for el, amt in composition.items():
            if el not in cs.list_of_elements|elements_variables.keys()|{"□"}:
                elements_variables[el] = []
            for var in re.findall("[a-z" + "".join(cnst.GREEK_CHARS) + "]", amt):
                stoichiometry_variables[var] = []

        rename_variables = [("R", "E"), ("A", "E"), ("T", "M")]
        for v1, v2 in rename_variables:
            if v1 in elements_variables and v2 in elements_variables and v1 + v2 in formula:
                elements_variables[v1 + v2] = []
                del elements_variables[v2]
                del elements_variables[v1]
                composition[v1 + v2] = composition[v2]
                del composition[v1]
                del composition[v2]

        if "M" in elements_variables and "e" in stoichiometry_variables:
            elements_variables["Me"] = []
            del elements_variables["M"]
            del stoichiometry_variables["e"]
            c = composition["M"][1:]
            composition["Me"] = c if c != "" else "1.0"
            del composition["M"]

        if not oxygen_deficiency and oxygen_deficiency_sym in stoichiometry_variables:
            oxygen_deficiency = None
        variables = [v for v in stoichiometry_variables.keys()
                     if [e for e, s in composition.items() if v in s] == ["O"]]
        oxygen_deficiency = variables[0] if len(variables) > 0 else oxygen_deficiency
        for var in variables:
            del stoichiometry_variables[var]
            composition["O"] = "1" if composition["O"] == var else composition["O"].replace(var, "").strip()
            formula = formula.replace(var, "")

        formula_structure = {"elements": composition,
                             "amounts_vars": {x: v for x, v in stoichiometry_variables.items()},
                             "elements_vars": {e: v for e, v in elements_variables.items()},
                             "formula": formula,
                             "oxygen_deficiency": oxygen_deficiency}

        return formula_structure

    def get_species(self, formula):
        number_to_alphabet_dict = {
            "specie0_": "A",
            "specie1_": "B",
            "specie2_": "C",
            "specie3_": "Q",
            "specie4_": "K",
            "specie5_": "F",
            "specie6_": "G",
            "specie7_": "H",
            "specie8_": "I",
            "specie9_": "J",
        }
        species_in_material, species_indexs, species_dict = {}, {}, collections.OrderedDict()
        material_formula = formula
        i = 0
        for species in cs.species:
            while species in material_formula:
                # print(species)
                material_formula = material_formula.replace(species, "specie" + str(i) + "_")
                species_in_material["specie" + str(i) + "_"] = species
                i += 1
        if species_in_material == {}:
            return None
        for species in number_to_alphabet_dict:
            while species in material_formula:
                material_formula = material_formula.replace(species, number_to_alphabet_dict[species])
                species_indexs[number_to_alphabet_dict[species]] = species_in_material[species]
        species_info = self.__parse_formula(material_formula)
        for species_index in species_info:
            species_dict[species_indexs[species_index]] = species_info[species_index]
        return species_dict

    def __fraction_convertion(self, formula):
        formula_upd = formula
        r_a = r"([0-9\.]*)"
        r_b = r"(\([0-9\.]*)"
        r_x = r"([a-z]*)"
        r_s = r"([\-\+]+)"
        r_d = r"([0-9\.]*)"
        r_y = r"([a-z]+\))"
        r_c = r"(?=[/]*([0-9\.]*))"
        r = r"(" + r_a + r_b + r_x + r_s + r_d + r_y + r_c + ")"
        for m in re.finditer(r, formula_upd):
            expr_old = m.group(1) + "/" + m.group(8) if m.group(8) != "" else m.group(1)
            a = m.group(2).strip(")(") if m.group(2).strip(")(") != '' else '1'
            b = m.group(3).strip(')(') if m.group(3).strip(')(') != '' else '1'
            x = m.group(4).strip(")(") if m.group(4).strip(")(") != '' else '1'
            s = m.group(5).strip(")(") if m.group(5).strip(")(") != '' else '+'
            d = m.group(6).strip(")(") if m.group(6).strip(")(") != '' else '1'
            y = m.group(7).strip(")(") if m.group(7).strip(")(") != '' else '1'
            c = m.group(8).strip(")(") if m.group(8).strip(")(") != '' else '1'
            expr_str = a + '/' + c + '*' + b + '*' + x + s + a + '/' + c + '*' + d + "*" + y
            expr = str(smp.simplify(expr_str)).strip()
            if expr[0] == '-':
                s_expr = re.split(r"\+", expr)
                expr = s_expr[1] + s_expr[0]
            expr_new = expr.strip().replace(" ", "")
            formula_upd = formula_upd.replace(expr_old, expr_new.strip(), 1)

        return re.sub(r"\s{1,}", "", formula_upd)

    def __get_oxygen_deficiency(self, formula):
        formula_upd = formula
        oxy_def = None
        oxy_def_sym = ""
        r = "".join([s for s in cnst.GREEK_CHARS])
        r = "O[0-9]*([-+±∓]{1})[a-z" + r + "]{1}[0-9]*$"
        for m in re.finditer(r, formula_upd.rstrip(")")):
            end = formula_upd[m.start():m.end()]
            splt = re.split("[-+±∓]", end)
            oxy_def_sym = splt[-1]
            oxy_def = m.group(1)
            formula_upd = formula_upd[:m.start()] + formula_upd[m.start():].replace(end, splt[0])

        if oxy_def_sym not in "[-+±∓]" and oxy_def_sym in re.findall("O(a-z)$", formula_upd.rstrip(")")):
            oxy_def_sym = "±"

        return formula_upd, oxy_def, oxy_def_sym

    def __parse_formula(self, init_formula):

        formula_dict = collections.OrderedDict()

        formula_dict = self.__parse_parentheses(init_formula, "1", formula_dict)
        """
        refinement of non-variable values
        """
        incorrect = []
        for el, amt in formula_dict.items():
            formula_dict[el] = self.__simplify(amt)
            if any(len(c) > 1 for c in re.findall("[A-Za-z]+", formula_dict[el])):
                incorrect.append(el)

        for el in incorrect:
            del formula_dict[el]

        return formula_dict

    def __parse_parentheses(self, init_formula, init_factor, curr_dict):
        r = r"\(((?>[^\(\)]+|(?R))*)\)\s*([-*\.\da-z\+/]*)"

        for m in re.finditer(r, init_formula):
            #print("--->", m.group(0), m.group(1), m.group(2))
            factor = "1"
            if m.group(2) != "":
                factor = m.group(2)

            factor = self.__simplify("(" + str(init_factor) + ")*(" + str(factor) + ")")
            unit_sym_dict = self.__parse_parentheses(m.group(1), factor, curr_dict)
            init_formula = init_formula.replace(m.group(0), "")

        unit_sym_dict = self.__get_sym_dict(init_formula, init_factor)
        for el, amt in unit_sym_dict.items():
            if el in curr_dict:
                if len(curr_dict[el]) == 0:
                    curr_dict[el] = amt
                else:
                    curr_dict[el] = "(" + str(curr_dict[el]) + ")" + "+" + "(" + str(amt) + ")"
            else:
                curr_dict[el] = amt

        return curr_dict

    def __get_sym_dict(self, f, factor):
        sym_dict = collections.OrderedDict()
        r = r"([A-Z□]{1}[a-z]{0,1})\s*([-\*\.\da-z" + "".join(cnst.GREEK_CHARS) + r"\+\/]*)"

        def get_code_value(code, iterator):
            code_mapping = {"01": (iterator.group(1), iterator.group(2)),
                            "11": (iterator.group(1), iterator.group(2)),
                            "10": (iterator.group(1)[0], iterator.group(1)[1:] + iterator.group(2)),
                            "00": (iterator.group(1)[0], iterator.group(1)[1:] + iterator.group(2))}
            return code_mapping[code]

        for m in re.finditer(r, f):
            """
            checking for correct elements names
            """
            el_bin = "{0}{1}".format(str(int(m.group(1)[0] in cs.list_of_elements_1|{"M", "□"})),
                                     str(int(m.group(1) in cs.list_of_elements|{"Ln", "M", "□"})))
            el, amt = get_code_value(el_bin, m)
            if amt.strip() == "":
                amt = "1"
            if el in sym_dict:
                sym_dict[el] = "(" + sym_dict[el] + ")" + "+" + "(" + amt + ")" + "*" + "(" + str(factor) + ")"
            else:
                sym_dict[el] = "(" + amt + ")" + "*" + "(" + str(factor) + ")"
            f = f.replace(m.group(), "", 1)
        if f.strip():
            return collections.OrderedDict()

        """
        refinement of non-variable values
        """
        try:
            for el, amt in sym_dict.items():
                sym_dict[el] = self.__simplify(amt)
        except:
            sym_dict = collections.OrderedDict()

        return sym_dict

    def __element_structure(self, element):
        return {"material_string": element,
                "material_name": cs.element2name.get(element.rstrip("0987654321"), ""),
                "material_formula": element,
                "additives": [],
                "phase": "",
                "oxygen_deficiency": None,
                "is_acronym": False,
                "amounts_vars": {},
                "elements_vars": {},
                "composition": [{"formula": element,
                                 "amount": "1",
                                 "elements": collections.OrderedDict([(element, "1")]),
                                 "species": collections.OrderedDict([(element, "1")])
                                 }]
                }

    ###################################################################################################################
    # Splitting list of materials
    ###################################################################################################################

    def is_materials_list(self, material_string):
        if (any(a + "s" in material_string.lower() for a in cs.anions.keys()) or "metal" in material_string) and \
                any(w in material_string for w in ["and ", ",", " of "]):
            return True

        return False

    def split_materials_list(self, material_string):
        """
        split material string into list of compounds when it"s given in form cation + several anions
        :param material_string: <str>
        :return: <list> of <str> chemical names
        """

        parts = [p for p in re.split(r"[\s\,]", material_string) if p != ""]

        anion = [(i, p[:-1]) for i, p in enumerate(parts)
                 if p[:-1].lower() in cs.anions.keys() or p[:-1].lower() == "metal"]
        cation = [(i, p) for i, p in enumerate(parts)
                  if p.lower() in cs.cations.keys() or p in cs.list_of_elements]
        valencies = [(i - 1, p.strip("()")) for i, p in enumerate(parts) if p.strip("()") in cs.rome2num and i != 0]

        result = []
        if len(anion) == 1:
            for c_i, c in cation:

                if c in cs.element2name:
                    name = [cs.element2name[c]]
                else:
                    name = [c.lower()]
                valency = "".join([v for v_i, v in valencies if v_i == c_i])
                if valency != "":
                    name.append("(" + valency + ")")
                name.append(anion[0][1])
                hydr_i = material_string.find("hydrate")
                if hydr_i > -1:
                    pref = []
                    while material_string[hydr_i - 1] != " " and hydr_i > 0:
                        pref.append(material_string[hydr_i - 1])
                        hydr_i -= 1

                    pref = "".join([p for p in reversed(pref)])

                    if pref not in cs.neg_prefixes:
                        name.append(pref + "hydrate")
                # result.append(("".join([n + " " for n in name]).strip(" "), valency))
                result.append(" ".join([n for n in name]))

        return result

    ###################################################################################################################
    # Reconstruct formula / dictionary lookup
    ###################################################################################################################

    def string2formula(self, material_string):

        if self.__verbose:
            print("Converting string to formula:")

        if material_string == '':
            return "", "", ""

        material_string_nodashes = ' '.join([t for t in re.split(r'\s|(?<=[a-z])-', material_string)])
        material_string_nodashes = material_string_nodashes[0].lower() + material_string_nodashes[1:]
        if material_string_nodashes in cs.pubchem_dictionary:
            return material_string, cs.pubchem_dictionary[material_string_nodashes], material_string_nodashes

        if material_string in cs.default_abbreviations:
            if self.__verbose:
                print("\tFound abbreviation:", material_string, "-->", cs.default_abbreviations[material_string])
            return material_string, cs.default_abbreviations[material_string], ""

        material_formula, material_name = self.extract_formula_from_string(material_string)
        if self.__verbose:
            print("\tAfter material name parsing into formula|name:\n",
                  material_string, "-->", material_formula, "|", material_name)

        if material_formula == "":
            material_formula = self.reconstruct_formula_from_string(material_name)
            if material_formula == "" and self.__pubchem:
                pcp_compounds = pcp.get_compounds(material_name, 'name')
                material_formula = pcp_compounds[0].molecular_formula if len(pcp_compounds) != 0 else ""
        if self.__verbose:
            print("\tAfter material formula reconstruction:\n",
                  material_string, "-->", material_formula, "|", material_name)

        material_name = "" if material_formula == "" else material_name
        material_formula = material_string if material_formula == "" else material_formula

        return material_string, material_formula, material_name

    def extract_formula_from_string(self, material_string):
        """
        Extracts chemical formula among chemical terms.
        This happens due to poor tokenization
        :param material_string
        :return: formula: <str> chemical formula found in material string
                 name: <str> chemical name found in material string
        """

        material_string = material_string.replace("[", " [").strip()
        material_string = re.sub(r"\s{2,}", " ", material_string)
        split = re.split(r"\s", material_string)

        if len(split) == 1:
            return "", material_string

        if len(split) == 2 and \
                all(t in cs.list_of_elements or t.rstrip('s').lower() in cs.list_of_anions|{'metal'}
                    for t in split):
            t1 = cs.element2name[split[0]] if split[0] in cs.element2name else split[0].rstrip('s').lower()
            t2 = cs.element2name[split[1]] if split[1] in cs.element2name else split[1].rstrip('s').lower()
            return "", t1 + " " + t2

        def try_formula(self, formula):
            if re.match(r"(\s*\([IV,]+\))", formula):
                return {}
            try:
                composition = self.formula2composition(formula)
            except:
                composition = {'elements': {}, 'elements_vars': {}}
            if self.__is_abbreviation_composition(composition):
                return {}
            return composition["elements"]

        formulas = []
        terms = []
        for part in split:
            composition_part = try_formula(self, part)
            if composition_part != {} or part.strip('+-1234567890') in cs.list_of_elements:
                formulas.append(part)
            else:
                terms.append(part)

        if len(formulas) > 1:
            # print("More than one formula")
            return material_string, ""

        if len(formulas) == 0:
            # print("No formula")
            return "", material_string
            # return "",  " ".join([t + " " for t in terms])

        if any(t.strip('+-1234567890') in cs.list_of_elements for t in terms):
            # print("Potentially many pieces of formula")
            return material_string, ""

        # print("One formula and chemical terms", terms)
        return formulas[0], " ".join([t for t in terms])

    def reconstruct_formula_from_string(self, material_name, valency=""):
        """
        reconstructing chemical formula for simple chemical names anion + cation
        :param material_name: <str> chemical name
        :param valency: <str> anion valency
        :return: <str> chemical formula
        """

        output_formula = ""
        material_name = re.sub(r"(\([IV]*\))", " \\1 ", material_name)
        material_name = re.sub(r"\s{2,}", " ", material_name)

        terms_list = []
        valency_list = []
        hydrate = ""
        cation_prefix_num = 0
        cation_data = {"c_name": "", "valency": [], "e_name": "", "n_atoms": 0}

        for t in material_name.split(" "):
            if t.strip("()") in cs.rome2num:
                valency_list.append(cs.rome2num[t.strip("()")])
                continue
            if "hydrate" in t.lower():
                hydrate = t
                continue
            terms_list.append(t.strip(" -"))

        terms_list_upd = []
        for t in terms_list:
            if all(p not in cs.prefixes2num for p in t.split("-")):
                terms_list_upd.extend(_ for _ in t.split("-"))
            else:
                terms_list_upd.append("".join([p for p in t.split("-")]))
        terms_list = terms_list_upd

        t = "".join([t + " " for t in terms_list]).lower().strip(" ")
        if t in cs.anions:
            return cs.anions[t]["e_name"]
        if t in cs.cations:
            return cs.cations[t]["e_name"]

        if len(terms_list) < 2:
            return output_formula

        anion = terms_list.pop().lower().rstrip("s")

        if valency == "":
            valency_num = max(valency_list + [0])
        else:
            valency_num = cs.rome2num[valency.strip("()")]

        next_term = terms_list.pop()
        if "hydrogen" in next_term.lower() and len(terms_list) != 0:
            anion = next_term + " " + anion
        else:
            terms_list += [next_term]

        _, anion_prefix_num, anion = ("", 0, anion) if anion.lower() in cs.anions else self.__get_prefix(anion)

        if anion in cs.anions:
            anion_data = cs.anions[anion].copy()
        elif anion in ["metal"] and terms_list[0].lower() in cs.cations:
            return cs.cations[terms_list[0].lower()]["e_name"]
        else:
            return output_formula

        if len(terms_list) >= 2:
            return output_formula

        if len(terms_list) == 1:
            cation = terms_list[0]
            _, cation_prefix_num, cation = ("", 0, cation) if cation.lower() in cs.cations \
                else self.__get_prefix(cation)
            if cation.lower() in cs.cations:
                cation = cation.lower()
                cation_data = cs.cations[cation].copy()
            elif cation in cs.element2name:
                cation_data = cs.cations[cs.element2name[cation]].copy()
            else:
                return output_formula

        if len(cation_data["valency"]) > 1 and valency_num != 0:
            if valency_num not in cation_data["valency"]:
                print("WARNING! Not common valency value for " + material_name)
                print(cation_data["valency"])
                print(valency_num)
            cation_data["valency"] = [valency_num]

        output_formula = self.__build_formula(anion=anion_data,
                                              cation=cation_data,
                                              cation_prefix_num=cation_prefix_num,
                                              anion_prefix_num=anion_prefix_num)

        if hydrate != "":
            _, hydrate_prefix_num, hydrate = self.__get_prefix(hydrate)
            hydrate_prefix = "" if hydrate_prefix_num in [0, 1] else str(hydrate_prefix_num)
            output_formula = output_formula + "·" + hydrate_prefix + "H2O"

        return output_formula

    def __build_formula(self, cation, anion, cation_prefix_num=0, anion_prefix_num=0):

        cation_stoich = 0
        anion_stoich = 0

        if anion_prefix_num + cation_prefix_num == 0 or anion_prefix_num * cation_prefix_num != 0:
            v_cation = abs(cation["valency"][0])
            v_anion = abs(anion["valency"][0])
            cm = self.__lcm(v_cation, v_anion)
            cation_stoich = cm // v_cation
            anion_stoich = cm // v_anion

        if anion_prefix_num != 0:
            anion_stoich = anion_prefix_num
            i = 0
            cation_stoich = 0
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

        return "{0}{1}{2}{3}".format(cation_name_el, self.__cast_stoichiometry(cation_stoich), anion_name_el,
                                     self.__cast_stoichiometry(anion_stoich))

    def __get_prefix(self, material_name):

        pref = ""
        pref_num = 0
        material_name_upd = material_name

        for p in cs.prefixes2num.keys():
            if material_name.lower().find(p) == 0 and p != "":
                pref = p
                pref_num = cs.prefixes2num[p]
                material_name_upd = material_name_upd[len(p):].strip("-")

                if material_name_upd == "xide":
                    material_name_upd = "oxide"

        return pref, pref_num, material_name_upd

    ###################################################################################################################
    # Splitting mixtures
    ###################################################################################################################

    def split_formula_into_compounds(self, material_name):
        """
        splitting mixture/composite/solid solution/alloy into compounds with fractions
        :param material_name: <str> material formula
        :return: <list> of <tuples>: (compound, fraction)
        """

        if material_name == "":
            return []

        if self.__verbose:
            print("Splitting formula into compounds:")

        split = self.__split_formula(material_name)
        l = 0
        while len(split) != l:
            l = len(split)
            split = [p for s in split for p in self.__split_formula(s[0], s[1])]

        output = []
        for m, f in split:
            try:
                f = smp.simplify(f)
                if f.is_Number:
                    f = round(float(f), 3)
                f = str(f)
            except:
                f = "1"

            f = self.__simplify(f)
            output.append((m, f))

        output = [(check_parentheses(m), f) for m, f in output]

        return output

    def __split_formula(self, material_name_, init_fraction="1"):

        re_str = r"(?<=[0-9\)])[-⋅·∙\∗](?=[\(0-9](?!x))|" + \
                 r"(?<=[A-Z])[-⋅·∙\∗](?=[\(0-9])|" + \
                 r"(?<=[A-Z\)])[-⋅·∙\∗](?=[A-Z])|" + \
                 r"(?<=[0-9\)])[-⋅·∙\∗](?=[A-Z])"
        re_str = re_str + "".join([r"|(?<=" + e + r")[-⋅·∙\∗](?=[\(0-9A-Z])" for e in cs.list_of_elements])
        re_str = re_str + r"|[-·]([nx0-9\.]H2O)"

        material_name = material_name_.replace(" ", "")
        # if "(1-x)" == material_name[0:5] or "(100-x)" == material_name[0:7]:
        #     material_name = material_name.replace("(x)", "x")
        #     parts = re.findall(r"\(10{0,2}-x\)(.*)[-+·∙\∗⋅]x(.*)", material_name)
        #     parts = parts[0] if parts != [] else (material_name[5:], "")
        #     return [(parts[0].lstrip(" ·*⋅"), "1-x"), (parts[1].lstrip(" ·*"), "x")]
        pref = [s for s in re.split("(^\(1-[xyz][-xyz]*\))|(^\(100-[xyz][-xyz]*\))", material_name) if s]
        if len(pref) > 1:
            material_name_temp = pref.pop()
            amount = pref.pop()
            variables = re.findall("[a-z]", amount)
            for v in variables:
                material_name = material_name.replace("(" + v + ")", v)
            compounds = []
            while variables:
                v = variables.pop()
                parts = re.findall(r"(.*)[-+·∙\∗⋅]" + v + "(.*)$", material_name_temp)
                if parts:
                    compounds.append((parts[0][1], v))
                    material_name_temp = parts[0][0]
            compounds.append((material_name_temp, amount.strip("()")))
            return [c for c in reversed(compounds)]

        parts = [p for p in re.split(re_str, material_name) if p]
        #print("-->", material_name, parts)

        if len(parts) > 1:
            parts_upd = [p for part in parts for p in
                         re.split(r"(?<=[A-Z\)])[-·∙\∗⋅](?=[xyz])|(?<=O[0-9\)]+)[-·∙\∗⋅](?=[xyz])", part)]
        else:
            parts_upd = parts

        if any(m.strip("0987654321") in cs.list_of_elements for m in parts_upd[:-1]):
            parts_upd = ["".join([p + "-" for p in parts_upd]).rstrip("-")]

        merged_parts = [parts_upd[0]]
        for m in parts_upd[1:]:
            if re.findall("[A-Z]", m) == ["O"]:
                to_merge = merged_parts.pop() + "-" + m
                merged_parts.append(to_merge)
            else:
                merged_parts.append(m)

        composition = []
        for m in merged_parts:
            fraction = ""
            i = 0
            while i < len(m) and not m[i].isupper():
                fraction = fraction + m[i]
                i += 1
            fraction = fraction.strip("()")
            if fraction == "":
                fraction = "1"
            else:
                m = m[i:]

            fraction = "(" + fraction + ")*(" + init_fraction + ")"

            if m != "":
                composition.append((m, fraction))
        return composition

    def separate_phase(self, material_name):
        phase = ""
        string_start = 0
        for m in re.finditer("^([A-Za-z" + "".join(cnst.GREEK_CHARS) + "][0-9]{0,1})\-[A-Z]\.*", material_name):
            phase = m.group(1)
            string_start = m.end()-1

        return phase, material_name[string_start:]

    def separate_additives(self, material_name):
        """
        resolving doped part in material string
        :param material_name: <str> material string
        :return: <list> of additives,
                <str> new material name
        """
        new_material_name = material_name
        additives = []

        new_material_name = new_material_name.replace("codoped", "doped")
        new_material_name = new_material_name.replace("co-doped", "doped")

        # checking for "doped with"
        for r in ("activated", "modified", "stabilized", "doped", "added"):
            parts = [w for w in re.split(r + " with", new_material_name) if w != ""]
            if len(parts) > 1:
                new_material_name = parts[0].strip(" -+")
                additives.append(parts[1].strip())

        # checking for element-doped prefix
        for r in ("activated", "modified", "stabilized", "doped", "added"):
            parts = [w for w in re.split(r"(.*)[-\s]{1}" + r + " (.*)", new_material_name) if w != ""]
            if len(parts) > 1:
                new_material_name = parts.pop()
                additives.extend(parts)

        if "%" in new_material_name:
            new_material_name = new_material_name.replace(".%", "%")
            parts = re.split(r"[\-+:·]{0,1}\s*[0-9x\.]*\s*[vmolwt\s]*\%", new_material_name)

        if len(parts) > 1:
            new_material_name = parts[0].strip(" -+")
            additives.extend(d.strip() for d in parts[1:] if d != "")

        for part_ in new_material_name.split(":"):
            part_ = part_.strip()

            part = part_
            if any(e in part for e in cs.list_of_elements_2):
                for e in cs.list_of_elements_2:
                    part = part.replace(e, "&&")

            if all(e.strip("zyx,+0987654321. ") in cs.list_of_elements_1|{"R", "&&"}
                   for e in re.split(r"[\s,/]", part) if e != ""):
                additives.append(part_.strip(" "))
            else:
                new_material_name = part_.strip(" ")

        additives_final = [a.strip() for s in additives for a in re.split(r"[\s,\-/]|and", s) if a.strip() != ""]

        return additives_final, new_material_name.rstrip("( ,.:;-±/+")

    def substitute_additives(self, additives, material_structure):
        """
        analyzes additives and adjust them to the composition
        :param additives: list of additives
        :param material_structure: parsed material structure
        :return: updated material structure
        """
        additive = additives[0]

        if self.__verbose:
            print("Substituting additives for:")
            pprint(material_structure)
            print(additives)

        material_structure_new = material_structure.copy()
        try:
            additive_composition = self.formula2composition(additive)
        except:
            additive_composition = self.__empty_composition().copy()
        if len(additive_composition["elements"]) > 1:
            # print('-->', "Additive is compound")
            for structure in material_structure_new["composition"]:
                structure["amount"] = structure["amount"] + "-x"
            material_structure_new["composition"].append(
                {"formula": additive,
                 "amount": "x",
                 "elements": additive_composition["elements"],
                 "species": self.get_species(additive)
                 }
            )
        elif all(c["elements"] != {} for c in material_structure_new["composition"]):
            # print('-->', "Additive is element with fraction")
            formula, composition = self.__substitute_additive(additive, material_structure_new["material_formula"],
                                                              material_structure_new["composition"])
            if formula != material_structure_new["material_formula"]:
                material_structure_new["additives"] = []
                material_structure_new["material_formula"] = formula
                material_structure_new["composition"] = composition
            else:
                material_structure_new["additives"] = additives
        else:
            # print('-->', "Additive is something else")
            material_structure_new["additives"] = additives

        return material_structure_new

    def __substitute_additive(self, additive, material_formula, material_composition):

        new_material_composition = []
        new_material_formula = material_formula

        if additive[-1] == "+":
            additive = additive.rstrip("+0987654321")

        r = r"^[x0-9\.]+|[x0-9\.]+$"
        coeff = re.findall(r, additive)
        element = [s for s in re.split(r, additive) if s != ""][0]

        if coeff == [] or element not in cs.list_of_elements:
            return new_material_formula, material_composition

        for compound in material_composition:
            expr = "".join(["(" + v + ")+" for e, v in compound["elements"].items()]).rstrip("+")

            coeff = coeff[0] if not re.match("^[0]+[1-9]", coeff[0]) else "0." + coeff[0][1:]
            expr = expr + "+(" + coeff + ")"

            if self.__is_int(self.__simplify(expr)):
                new_name = element + coeff + compound["formula"]
                new_composition = compound["elements"].copy()
                new_composition.update({element: coeff})
                new_composition.move_to_end(element, last=False)

                new_material_composition.append({"formula": new_name,
                                                 "amount": compound.get("amount"),
                                                 "elements": new_composition
                                                 })
                new_material_formula = new_material_formula.replace(compound.get("formula"), new_name)
            else:
                new_material_composition.append({"formula": compound.get("formula"),
                                                 "amount": compound.get("amount"),
                                                 "elements": compound.get("elements")
                                                 })

        return new_material_formula, new_material_composition

    ###################################################################################################################
    # Resolving abbreviations
    ###################################################################################################################

    def __is_abbreviation(self, word):
        if all(c.isupper() for c in re.sub(r"[0-9x\-\(\)\.]", "", word)) and len(re.findall("[A-NP-Z]", word)) > 1:
            return True

        return False

    def build_acronyms_dict(self, materials_list, paragraph):
        """
        constructing dictionary of acronyms appeared in material list
        :param paragraph: <list> of <str> list of sentences to look for abbreviations names
        :param materials_list: <list> of <str> list of materials entities
        :return: <dict> abbreviation: corresponding string
        """

        acronyms_dict = {t: "" for t in materials_list if self.__is_abbreviation(t.replace(" ", "")) and t != ""}
        not_abbreviations = list(set(materials_list) - set(acronyms_dict.keys()))

        # first find acronyms in current materials list
        for abbr in acronyms_dict.keys():
            for material_name in not_abbreviations:
                if sorted(re.findall("[A-NP-Z]", abbr)) == sorted(re.findall("[A-NP-Z]", material_name)):
                    acronyms_dict[abbr] = material_name

        # for all other acronyms going through the paper text
        for abbr, name in acronyms_dict.items():
            sents = " ".join([s + " " for s in paragraph if abbr in s]).strip(" ").split(abbr)
            i = 0
            while acronyms_dict[abbr] == "" and i < len(sents):
                sent = sents[i]
                for tok in sent.split(" "):
                    if sorted(re.findall("[A-NP-Z]", tok)) == sorted(re.findall("[A-NP-Z]", abbr)):
                        acronyms_dict[abbr] = tok
                i += 1

        for abbr in acronyms_dict.keys():
            parts = re.split("-", abbr)
            if all(p in acronyms_dict for p in parts) and acronyms_dict[abbr] == "" and len(parts) > 1:
                name = "".join("(" + acronyms_dict[p] + ")" + "-" for p in parts).rstrip("-")
                acronyms_dict[abbr] = name

        empty_list = [abbr for abbr, name in acronyms_dict.items() if name == ""]
        for abbr in empty_list:
            del acronyms_dict[abbr]

        return acronyms_dict

    ###################################################################################################################
    # Methods to substitute variables
    ###################################################################################################################

    def __get_values(self, string, mode):
        values = []
        max_value = None
        min_value = None

        if len(string) == 0:
            return {"values": [],
                    "max_value": None,
                    "min_value": None}

        # given range
        if mode == "range":
            min_value, max_value = string[0]
            max_value = max_value.rstrip("., ")
            min_value = min_value.rstrip("., ")
            max_value = re.sub("[a-z]*", "", max_value)
            min_value = re.sub("[a-z]*", "", min_value)
            try:
                max_value = round(float(smp.simplify(max_value)), 4)
                min_value = round(float(smp.simplify(min_value)), 4) if min_value != "" else None
            except Exception as ex:
                max_value = None
                min_value = None
                template = "An exception of type {0} occurred when use sympy. Arguments:\n{1!r}."
                message = template.format(type(ex).__name__, ex.args)
                print(message)
            values = []

            return {"values": values,
                    "max_value": max_value,
                    "min_value": min_value}

        # given list
        if mode == "values":
            values = re.split(r"[,\s]", re.sub("[a-z]+", "", string[0]))
            try:
                values = [round(float(smp.simplify(c.rstrip("., "))), 4) for c in values if
                          c.rstrip("., ") not in ["", "and"]]
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

        regs = [(var + r"\s*=\s*([-]{0,1}[0-9\.\,/and\s]+)[\s\)\]\,]", "values"),
                (var + r"\s*=\s*([0-9\.]+)\s*[-–]\s*([0-9\.\s]+)[\s\)\]\,m\%]", "range"),
                (r"([0-9\.\s]*)\s*[<≤⩽]{0,1}\s*" + var + r"\s*[<≤⩽>]{1}\s*([0-9\.\s]+)[\s\)\]\.\,]", "range"),
                (var + r"[a-z\s]*from\s([0-9\./]+)\sto\s([0-9\./]+)", "range"),
                ]

        for r, m in regs:
            if values["values"] == [] and values["max_value"] is None:
                r_res = re.findall(r, sentence.replace(" - ", "-"))
                values = self.__get_values(r_res, m)

        return values

    def get_elements_values(self, var, sentence):
        """
        find elements values for var in the sentence
        :param var: <str> variable name
        :param sentence: <str> sentence to look for
        :return: <list> of <str> found values
        """
        values = re.findall(var + r"\s*[=:]{1}\s*([A-Za-z0-9\+,\s]+)", sentence)
        values = [c.rstrip("0987654321+") for v in values for c in re.split(r"[,\s]", v)
                  if c.rstrip("0987654321+") in cs.list_of_elements]

        return list(set(values))

    def __get_substitutions_array(self, subs_dict):

        """
        Generates combinations of different variables values
        I.e. if 'x' = [0.1, 0.2] and 'y' = [0.5, 0.6], then outputs: [
        {'x': 0.1, 'y': 0.5}, {'x': 0.1, 'y': 0.6}, {'x': 0.2, 'y': 0.5},  {'x': 0.2, 'y': 0.6}]
        :param subs_dict: dict(var: list of values)
        :return: list of dict(var: value)
        """
        subs_array = []
        l_dict = len(subs_dict)
        t_array = [{"var": k, "val": v} for k, vs in subs_dict.items() for v in vs["values"]]
        for comb in itertools.combinations(range(0, len(t_array)), l_dict):
            s = ''.join([t_array[i]['var'] for i in comb])
            if len(s) == len(set(s)):
                t_dict = {}
                for i in comb:
                    t_dict[t_array[i]['var']] = t_array[i]['val']
                subs_array.append(t_dict)

        return subs_array

    def substitute_amounts(self, material_structure):
        """
        substituting values for elements fractions variables into formula
        :param material_structure: <dict> output of mp.parse_material() with filled "fraction_vars"
        :return: list of structures derived from input with substitution of all fraction_vars
        """
        def update_stoichiometry(stoich, substitution):
            for var, val in subs.items():
                stoich_upd = stoich.replace(var, str(val))
            try:
                stoich_upd = self.__simplify(stoich_upd)
                if stoich_upd[0] == "-" and re.findall("[a-z]", stoich_upd) == []:
                    stoich_upd = "NEG"
            except:
                stoich_upd = stoich
            return stoich_upd

        new_materials_array = []
        fraction_variables = {x: v['values'] for x, v in material_structure['amounts_vars'].items()}
        fractions_array = self.__get_substitutions_array(fraction_variables)

        for subs in fractions_array:
            material_composition = []
            for m in material_structure['composition']:
                composition = {"formula": m.get('formula'),
                               "amounts": m.get('amount'),
                               "elements": collections.OrderedDict(),
                               "species": collections.OrderedDict()}

                for el, stoich in m['elements'].items():
                    stoich_upd = update_stoichiometry(stoich, subs)
                    if stoich_upd not in ["0.0", "0"]:
                        composition['elements'][el] = stoich_upd

                material_composition.append(composition)

            if all(v != "NEG" for c in material_composition for e, v in c["elements"].items()):
                material_structure_upd = material_structure.copy()
                material_structure_upd['elements_vars'] = material_structure['elements_vars'].copy()
                material_structure_upd["composition"] = material_composition
                material_structure_upd["amount_vars"] = subs
                new_materials_array.append(material_structure_upd)

        return new_materials_array


    ###################################################################################################################
    # Misc
    ###################################################################################################################


    def __simplify(self, value):
        """
        simplifying stoichiometric expression
        :param value: string
        :return: string
        """
        for l in cnst.GREEK_CHARS:
            _clash[l] = smp.Symbol(l)

        new_value = value
        for i, m in enumerate(re.finditer(r"(?<=[0-9])([a-z" + "".join(cnst.GREEK_CHARS) + "])", new_value)):
            new_value = new_value[0:m.start(1) + i] + "*" + new_value[m.start(1) + i:]
        new_value = smp.simplify(smp.sympify(new_value, _clash))
        if new_value.is_Number:
            new_value = round(float(new_value), 3)
        else:
            new_value = new_value.evalf(3)
        new_value = re.sub('\.0+(?![0-9])', '', str(new_value).replace(" ", ""))
        if new_value[0] == "-":
            split = re.split("\+", new_value)
            new_value = split[1] + split[0] if len(split) == 2 else new_value

        return new_value

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

        value = float(value)
        if value == 1.0:
            return ""
        if value * 1000 % 1000 == 0.0:
            return str(int(value))

        return str(value)

    def __empty_composition(self):
        return {"formula": "",
                "elements": collections.OrderedDict(),
                "amounts_vars": {}, "elements_vars": {},
                "phase": None, "oxygen_deficiency": None}

    def __empty_structure(self):
        return {"material_string": "", "material_name": "", "material_formula": "",
                "phase": None, "additives": [], "oxygen_deficiency": None,
                "is_acronym": False,
                "amounts_vars": {}, "elements_vars": {},
                "composition": []}

    def __is_int(self, num):
        try:
            return round(float(num), 3) == round(float(num), 0)
        except:
            return False

    def __is_abbreviation_composition(self, composition):
        if all(e.isupper() and s in ["1.0", "1", "s"] for e, s in composition["elements"].items()):
            return True

        elements_vars = [el for el in composition["elements_vars"].keys() if len(el) == 1 and el.isupper()]
        if len(elements_vars) > 1:
            return True

        return False

    def __is_acronym(self, structure):
        if any(all(e.isupper() and s in ["1.0", "1"] for e, s in compound["elements"].items()) for compound in
               structure["composition"]):
            return True

        elements_vars = [el for el in structure["elements_vars"].keys() if len(el) == 1 and el.isupper()]
        if len(elements_vars) > 1:
            return True

        if all(c.isupper() for c in structure["material_formula"]) and any(
                c not in cs.list_of_elements_1 for c in structure["material_formula"]):
            return True

        if re.findall("[A-Z]{3,}", structure["material_formula"]) != [] and \
                all(w not in structure["material_formula"] for w in
                    ["CH", "COO", "OH", "NH"] + [a for a in cs.default_abbreviations.keys()]):
            return True

        if "PV" == structure["material_formula"][0:2]:
            return True

        return False

    def get_element(self, name):
        if name in cs.anions:
            return cs.anions[name]["e_name"]
        if name in cs.cations:
            return cs.cations[name]["e_name"]
        return ""

    def __combine_formula(self, material_composition):
        formula = ""
        if len(material_composition) == 1:
            return material_composition[0]["formula"].replace("*", "")

        coeff = ""
        for c in material_composition:
            if coeff != "":
                coeff = self.__simplify(coeff)
            if all(ch.isdigit() or ch == "." for ch in c["amount"]):
                coeff = self.__cast_stoichiometry(c["amount"])
            else:
                coeff = "(" + c["amount"] + ")" if len(c["amount"]) != 1  else c["amount"]

            sign = "-"
            if "H2O" in c["formula"]:
                sign = "·"

            formula = formula + sign + coeff + c["formula"]

        formula = formula.replace("*", "")

        return formula.lstrip("-")
