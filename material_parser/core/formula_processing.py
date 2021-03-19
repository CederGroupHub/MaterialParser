# coding=utf-8
from collections import defaultdict, OrderedDict
from material_parser.core.utils import simplify

import regex as re
import sympy as smp
import material_parser.core.regex_parser as rp
import material_parser.core.chemical_sets as cs
import material_parser.core.constants as cnst

from pprint import pprint


def separate_phase(formula):
    phase = ""
    start = 0
    for m in re.finditer(rp.re_phase_prefix, formula):
        phase = m.group(1)
        start = m.end() - 1

    return phase, formula[start:]


def process_formula(formula):
    """
    :param formula: str
    :return: chemical_structure
    """
    """
    build dictionary of elements and stoichiometries
    """
    formula_data = parse_formula(formula)

    """
    looking for variables in elements and stoichiometry
    """
    for el, amt in formula_data["composition"].items():
        if el not in cs.list_of_elements | formula_data["elements_x"].keys() | cnst.VACANCIES:
            formula_data["elements_x"][el] = []
        for var in re.findall(rp.re_variables, amt):
            formula_data["amounts_x"][var] = []

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

    if __is_acronym(formula, elements, elements_x):
        return dict(formula=formula,
                    elements=OrderedDict(),
                    species=OrderedDict(),
                    oxygen_deficiency="",
                    phase="",
                    amounts_x={},
                    elements_x={})

    species = __get_species(formula) if len(elements) > 2 or formula == "H2O" else elements

    return dict(formula=formula,
                elements=elements,
                species=species,
                oxygen_deficiency=oxygen_deficiency,
                phase=formula_data["phase"],
                amounts_x={x: v for x, v in stoichiometry_x.items()},
                elements_x={e: v for e, v in elements_x.items()})


def parse_formula(formula):
    formula = formula.replace(" ", "")
    """
    separate phase, e.g. g-ABC
    """
    phase, formula = separate_phase(formula)

    """
    separate oxygen deficiency
    """
    formula, oxygen_deficiency, oxygen_deficiency_sym = __get_oxygen_deficiency(formula)

    """
    converting fractions a(b+x)/c into (a/c*b+a/c*x)
    """
    formula = __fraction_convertion(formula)

    elements_x = defaultdict(str)
    stoichiometry_x = defaultdict(str)

    """
    check for any weird syntax (A,B)zElxEly...
    replacing with MzElxEly... and M = [A, B]
    """
    for m in re.finditer(rp.re_weird_syntax, formula):
        elements_x["M"] = re.split(r"[\/,]", m.group(0).strip('()'))
        formula = formula.replace(m.group(0), "M", 1)

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

    for m in re.finditer(rp.re_in_parentheses, init_formula):
        # print("--->", m.group(0), m.group(1), m.group(2))
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
    sym_dict = OrderedDict()

    def get_code_value(code, iterator):
        code_mapping = {"01": (iterator.group(1), iterator.group(2)),
                        "11": (iterator.group(1), iterator.group(2)),
                        "10": (iterator.group(1)[0], iterator.group(1)[1:] + iterator.group(2)),
                        "00": (iterator.group(1)[0], iterator.group(1)[1:] + iterator.group(2))}
        return code_mapping[code]

    for m in re.finditer(rp.re_sym_dict, f):
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


def __get_species(formula):
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
    species_info = parse_formula(material_formula)["composition"]
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


def __get_oxygen_deficiency(formula):
    formula_upd = formula
    oxy_def = ""
    oxy_def_sym = ""

    if len(formula) < 3:
        return formula_upd, oxy_def, oxy_def_sym

    if formula[-2:] in cs.list_of_elements_2:
        return formula_upd, oxy_def, oxy_def_sym

    for m in re.finditer(rp.re_oxygen_deficiency, formula_upd.rstrip(")")):
        end = formula_upd[m.start():m.end()]
        splt = re.split(rp.re_signs, end)
        oxy_def_sym = splt[-1]
        oxy_def = m.group(1)
        formula_upd = formula_upd[:m.start()] + formula_upd[m.start():].replace(end, splt[0])

    if oxy_def_sym not in rp.re_signs and oxy_def_sym == formula_upd.rstrip(")")[-1]:
        oxy_def_sym = "±"

    return formula_upd, oxy_def, oxy_def_sym


def __fraction_convertion(formula):
    formula_upd = formula
    for m in re.finditer(rp.re_formula_fraction, formula_upd):
        expr_old = m.group(1) + "/" + m.group(8) if m.group(8) != "" else m.group(1)
        a = m.group(2).strip(")(") if m.group(2).strip(")(") != '' else '1'
        b = m.group(3).strip(")(") if m.group(3).strip(')(') != '' else '1'
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


def __is_acronym(formula, composition, variables):
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
