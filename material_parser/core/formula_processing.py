# coding=utf-8
from collections import defaultdict, OrderedDict
from material_parser.core.utils import simplify

import regex as re
import sympy as smp
import material_parser.core.regex_parser as rp
import material_parser.core.chemical_sets as cs
import material_parser.core.constants as cnst

from pprint import pprint


def process_formula(formula, regex_parser):
    """
    :param regex_parser:
    :param formula: str
    :return: chemical_structure
    """
    for s, r in cs.species_acronyms.items():
        formula = formula.replace(s, r)
    """
    build dictionary of elements and stoichiometries
    """
    formula_data = parse_formula(formula, regex_parser)

    """
    looking for variables in elements and stoichiometry
    """
    for el, amt in formula_data["composition"].items():
        if el not in cs.list_of_elements | formula_data["elements_x"].keys() | cnst.VACANCIES:
            formula_data["elements_x"][el] = []
        for var in re.findall(rp.re_variables, amt):
            formula_data["amounts_x"][var] = {}

    formula, \
    elements, \
    elements_x, \
    stoichiometry_x, \
    oxygen_deficiency = __refine_variables(formula_data["formula"],
                                           formula_data["composition"],
                                           formula_data["elements_x"],
                                           formula_data["amounts_x"],
                                           formula_data["oxygen_deficiency"],
                                           formula_data["oxygen_deficiency_sym"])

    if __is_acronym(formula, elements, elements_x) or __has_negative_composition(elements):
        return dict(formula=formula,
                    elements=OrderedDict(),
                    species=OrderedDict(),
                    oxygen_deficiency="",
                    phase="",
                    amounts_x={},
                    elements_x={})

    species = __get_species(formula, regex_parser) if len(elements) > 2 or formula == "H2O" else elements

    # print (formula)
    # pprint(elements)

    return dict(formula=formula,
                elements=elements,
                species=species,
                oxygen_deficiency=oxygen_deficiency,
                phase=formula_data["phase"],
                amounts_x={x: v for x, v in stoichiometry_x.items()},
                elements_x={e: v for e, v in elements_x.items()})


def parse_formula(formula, regex_parser):
    formula = formula.replace(" ", "")
    """
    separate phase, e.g. g-ABC
    """
    phase, formula = regex_parser.separate_phase(formula)

    """
    separate oxygen deficiency
    """
    formula, oxygen_deficiency, oxygen_deficiency_sym = regex_parser.separate_oxygen_deficiency(formula)

    """
    converting fractions a(b+x)/c into (a/c*b+a/c*x)
    """
    formula = regex_parser.make_fraction_convertion(formula)

    """
    check for any weird syntax (A,B)zElxEly...
    replacing with MzElxEly... and M = [A, B]
    """
    elements_x = defaultdict(str)
    stoichiometry_x = defaultdict(str)
    formula, variables = regex_parser.convert_weird_syntax(formula)
    if variables:
        elements_x["M"] = variables

    composition = __get_composition(formula)

    return dict(formula=formula,
                composition=composition,
                oxygen_deficiency=oxygen_deficiency,
                oxygen_deficiency_sym = oxygen_deficiency_sym,
                phase=phase,
                amounts_x={x: v for x, v in stoichiometry_x.items()},
                elements_x={e: v for e, v in elements_x.items()})


def __get_composition(init_formula):
    """

    :param init_formula:
    :return:
    """
    """
    if more than 4 repeating lowercase letters encountered then it is not chemical formula
    """
    if re.findall("[a-z]{4,}", init_formula):
        return OrderedDict()

    formula_dict = OrderedDict()
    formula_dict = __parse_parentheses(init_formula, "1", formula_dict)

    """
    refinement of non-variable values
    """
    incorrect = []
    for el, amt in formula_dict.items():
        formula_dict[el] = simplify(amt)
        if any(len(c) > 1 for c in re.findall("[A-Za-z]+", formula_dict[el])):
            incorrect.append(el)

    for el in incorrect:
        del formula_dict[el]

    return formula_dict


def __parse_parentheses(init_formula, init_factor, curr_dict):
    re_in_parentheses = r"\(((?>[^\(\)]+|(?R))*)\)\s*([-*\.\da-z\+/]*)"
    for m in re.finditer(re_in_parentheses, init_formula):
        factor = m.group(2) if m.group(2) != "" else "1"
        factor = simplify("(" + str(init_factor) + ")*(" + str(factor) + ")")
        unit_sym_dict = __parse_parentheses(m.group(1), factor, curr_dict)
        init_formula = init_formula.replace(m.group(0), "")

    unit_sym_dict = __get_sym_dict(init_formula, init_factor)
    for el, amt in unit_sym_dict.items():
        if el in curr_dict:
            if len(curr_dict[el]) != 0:
                curr_dict[el] = "(" + str(curr_dict[el]) + ")" + "+" + "(" + str(amt) + ")"
            else:
                curr_dict[el] = amt
        else:
            curr_dict[el] = amt

    return curr_dict


def __get_sym_dict(f, factor):
    re_sym_dict = r"([A-Z□]{1}[a-z]{0,1})\s*([\-\*\.\da-z" + "".join(cnst.GREEK_CHARS) + r"\+\/]*)"
    sym_dict = OrderedDict()

    def get_code_value(code, iterator):
        code_mapping = {"01": (iterator.group(1), iterator.group(2)),
                        "11": (iterator.group(1), iterator.group(2)),
                        "10": (iterator.group(1)[0], iterator.group(1)[1:] + iterator.group(2)),
                        "00": (iterator.group(1)[0], iterator.group(1)[1:] + iterator.group(2))}
        return code_mapping[code]

    for m in re.finditer(re_sym_dict, f):
        """
        checking for correct elements names
        """
        el_bin = "{0}{1}".format(str(int(m.group(1)[0] in cs.list_of_elements_1 | {"M"} | cnst.VACANCIES)),
                                 str(int(m.group(1) in cs.list_of_elements | {"Ln", "M"} | cnst.VACANCIES)))
        el, amt = get_code_value(el_bin, m)
        if amt.strip() == "":
            amt = "1"
        if el in sym_dict:
            sym_dict[el] = "(" + sym_dict[el] + ")" + "+" + "(" + amt + ")" + "*" + "(" + str(factor) + ")"
        else:
            sym_dict[el] = "(" + amt + ")" + "*" + "(" + str(factor) + ")"
        f = f.replace(m.group(), "", 1)
    if f.strip():
        return OrderedDict()

    """
    refinement of non-variable values
    """
    try:
        for el, amt in sym_dict.items():
            sym_dict[el] = simplify(amt)
    except:
        sym_dict = OrderedDict()
    return sym_dict


def __get_species(formula, regex_parser):
    species_in_material, species_indexs, species_dict = OrderedDict(), OrderedDict(), OrderedDict()
    material_formula = formula
    i = 0
    for species in cs.species:
        while species in material_formula:
            # print(species)
            material_formula = material_formula.replace(species, "specie" + str(i) + "_")
            species_in_material["specie" + str(i) + "_"] = species
            i += 1

    if not species_in_material:
        return OrderedDict()

    for species in cs.number_to_alphabet_dict:
        while species in material_formula:
            material_formula = material_formula.replace(species, cs.number_to_alphabet_dict[species])
            species_indexs[cs.number_to_alphabet_dict[species]] = species_in_material[species]
    species_info = parse_formula(material_formula, regex_parser)["composition"]
    for species_index in species_info:
        species_dict[species_indexs[species_index]] = species_info[species_index]
    return species_dict


def __refine_variables(formula, composition, elements_vars, stoichiometry_vars, oxy_def, oxy_def_sym):
    """
    :return:
    """
    """
    combining [RE, AE, TM] into one variable
    """
    rename_variables = [("R", "E"), ("A", "E"), ("T", "M")]
    for v1, v2 in rename_variables:
        if v1 in elements_vars and v2 in elements_vars and v1 + v2 in formula:
            elements_vars[v1 + v2] = []
            del elements_vars[v2]
            del elements_vars[v1]
            composition[v1 + v2] = composition[v2]
            del composition[v1]
            del composition[v2]

    """
    correction for Me variable
    """
    if "M" in elements_vars and "e" in stoichiometry_vars:
        elements_vars["Me"] = []
        del elements_vars["M"]
        del stoichiometry_vars["e"]
        c = composition["M"][1:]
        composition["Me"] = c if c != "" else "1.0"
        del composition["M"]

    """
    remove oxygen deficiency from variables
    """
    if not oxy_def and oxy_def_sym in stoichiometry_vars:
        oxy_def = None
    variables = [v for v in stoichiometry_vars.keys()
                 if [e for e, s in composition.items() if v in s] == ["O"]]
    oxy_def = chr(177) if len(variables) > 0 else oxy_def
    for var in variables:
        del stoichiometry_vars[var]
        composition["O"] = "1" if composition["O"] == var else composition["O"].replace(var, "").strip()
        formula = formula.replace(var, "")
    return formula, composition, elements_vars, stoichiometry_vars, oxy_def


def __is_acronym(formula, composition, variables):

    if formula in cs.ions:
        return False

    if any(ion in formula and len(ion) > 1 for ion in cs.ions):
        return False

    if len(composition) == 2 and variables:
        return True

    capital_letters = cnst.LATIN_CAPITAL - set(cs.list_of_elements_1) - {"M", "L"}
    if [r for c in capital_letters for r in re.findall(c + "[A-Z0-9\-]", formula)] \
            and all(w not in formula for w in ["RE", "OAC", "TM", "ME"]):
        return True

    if all(e.isupper() and s in ["1.0", "1"] for e, s in composition.items()):
        return True

    elements_x = [el for el in variables.keys() if len(el) == 1 and el.isupper()]
    if len(elements_x) > 1:
        return True

    if all(c.isupper() for c in formula) and any(c not in cs.list_of_elements_1 for c in formula):
        return True

    if re.findall("[A-Z]{3,}", formula) != [] and \
        all(w not in formula for w in ["CH", "COO", "OH", "NH"] + [a for a in cs.default_abbreviations.keys()]):
        return True

    if "PV" == formula[0:2]:
        return True

    return False


def __has_negative_composition(composition):

    flag = False
    try:
        flag = any(float(amt) < 0 for el, amt in composition.items())
    except:
        pass

    return flag
