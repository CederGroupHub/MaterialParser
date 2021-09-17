# coding=utf-8
import json
import os
import regex as re
import material_parser.core.chemical_sets as cs
import material_parser.core.constants as C
import material_parser.core.regex_parser as rp

from material_parser.core.utils import check_parentheses, parentheses_balanced


RE_COMPARE_SIGNS = r"\({0,1}[0-9\.]*\s*[" + \
                   "".join(C.COMPARE_SIGNS) + \
                   r"]{0,1}\s*[x,y]{0,1}\s*[" + \
                   "".join(C.COMPARE_SIGNS) + \
                   r"=]\s*[0-9\.-]*\){0,1}"

RE_OXY_DEF_REPLACE = r"O[0-9\.]*\s*[±\+\-∓]\s*([" + r"".join(C.GREEK_CHARS) + r"]{1})"

RE_TAIL_TRASH = [r"\(⩾99", r"\(99", r"\(98", r"\(90", r"\(95", r"\(96", r"\(Alfa", r"\(Aldrich", r"\(A.R.",
                 r"\(Aladdin", r"\(Sigma", r"\(A.G", r"\(Fuchen", r"\(Furuuchi", r"\(AR\)", "（x", r"\(x",
                 r"\(Acr[a-z]*", r"\(Koj", r"\(Sho", r"\(＞99"]

TRASH_SYMBOLS = ["#", "$", "!", "@", "©", "®", chr(8201), "Ⓡ", "\u200b"]

TRASH_TERMS = ["powder", "ceramic", "rear", "earth", "micro", "nano", "coat", "crystal", "particl", "glass"]

re_parentheses_open = ["{", "["]
re_parentheses_close = ["}", "]"]

__filename = os.path.dirname(os.path.realpath(__file__))
typos = json.loads(open(os.path.join(__filename, "rsc/typos.json")).read())


def combine_formula_parts(formula):
    """
    removes extra spaces if they break formula
    :param formula:
    :return:
    """
    parts = [p for p in reversed(re.split(r"\s", formula))]
    formula_upd = ""
    while parts:
        part = parts.pop()
        if parentheses_balanced(part):
            formula_upd = formula_upd + part + " "
        else:
            part_i = part
            while parts and not parentheses_balanced(part_i):
                part = parts.pop()
                part_i = part_i + part
            formula_upd = formula_upd + part_i + " "

    formula_upd = formula if formula_upd == "" else formula_upd.strip()

    return formula_upd


def cleanup_name(material_name):
    """
    cleaning up material name - fix due to tokenization imperfectness - TO BE REMOVED FROM FINAL VERSION
    :param material_name: <str> material string
    :return: <str> updated material string
    """
    # TODO: [Ti(N3)6]2−, Cu(NO3)2⋅4 H2O

    # correct dashes
    re_str = "\s*[" + "".join(C.DASHES) + "]\s*"
    material_name = re.sub(re_str, chr(45), material_name)
    material_name = re.sub(r"\s*\+\s*", "+", material_name)

    material_name = material_name.replace(chr(160), "")

    # correcting dots
    re_str = "[\\" + "".join(C.DOTS) + "]"
    material_name = re.sub(re_str, chr(183), material_name)
    material_name = re.sub(r"\s" + chr(183) + r"\s", chr(183), material_name)
    if chr(183) not in material_name:
        material_name = re.sub(r"\.([0-9x]*H2O)", chr(183)+"\\1", material_name)

    # correcting slashes
    slashes = [8725]
    re_str = "".join([chr(c) for c in slashes])
    re_str = "[\\" + re_str + "]"
    material_name = re.sub(re_str, chr(47), material_name)

    for c in re.findall(RE_OXY_DEF_REPLACE, material_name):
        material_name = material_name.replace(c, chr(948))

    # removing phase
    for c in ["(s)", "(l)", "(g)", "(aq)"]:
        material_name = material_name.replace(c, "")
    material_name = re.sub("\([a-w]{2}\.{0,1}\)$", "", material_name)

    # removing numbers in front of the formula
    material_name = re.sub("^[0-9]{3,}", "", material_name)

    # removing trash words
    for word in TRASH_TERMS:
        material_name = re.sub("[A-Za-z-]*" + word + "[a-z-]*", "", material_name)
        material_name = re.sub(word.capitalize() + "[a-z-]*", "", material_name)

    if any(a in material_name for a in C.ARROWS):
        return ""

    if "hbox" in material_name.lower():
        material_name = re.sub(r"(\\\\[a-z\(\)]+)", "", material_name)
        for t in ["{", "}", "_", " "]:
            material_name = material_name.replace(t, "")
        material_name = material_name.rstrip("\\")

    material_name = re.sub(RE_COMPARE_SIGNS, "", material_name)

    if material_name == "" or len([c for c in material_name if c.isalpha()]) < 1:
        return ""

    for c in RE_TAIL_TRASH:
        split = re.split(c, material_name)
        if len(split) > 1 and (len(split[-1]) == "" or all(not s.isalpha() for s in split[-1])):
            material_name = "".join([s for s in split[:-1]])

    material_name = re.sub("[C]{0,1}N[TFP]{1}[s]{0,1}", "", material_name)
    if material_name[-2:] not in ["Cs", "As", "Os"]:
        material_name = re.sub("[\-]([A-NP-Z]+[s]{0,1})$", "", material_name)

    for typo, correct in typos.items():
        material_name = material_name.replace(typo, correct)

    if material_name[-2:] == "/C":
        material_name = material_name[:-2]

    material_name = re.sub(r"(poly[\s-])(?=[a-z])", "poly", material_name)

    if material_name[-2:] == "/C":
        material_name = material_name[:-2]

    if len(material_name.split(" ")) > 1:
        for v in re.findall(r"[a-z](\([IV,]+\))", material_name):
            material_name = material_name.replace(v, " " + v)

    # if material_name != "":
    #     material_name = cs.check_parentheses(material_name)

    for c in TRASH_SYMBOLS:
        material_name = material_name.replace(c, "")

    for c in re_parentheses_open:
        material_name = material_name.replace(c, "(")
    for c in re_parentheses_close:
        material_name = material_name.replace(c, ")")

    material_name = material_name.lstrip(") -")
    material_name = material_name.rstrip("( ,.:;-±/∓")

    if len(material_name) == 1 and material_name not in cs.list_of_elements_1:
        return ""

    if len(material_name) == 2 and \
                    material_name not in cs.list_of_elements_2 and \
                    material_name.rstrip("234") not in cs.list_of_elements_1 and \
            any(c not in cs.list_of_elements_1 for c in material_name):
        return ""

    material_name = re.sub(r"\s([0-9\.]*H2O)$", chr(183) + "\\1", material_name)

    # print("-->", material_name)
    material_name = combine_formula_parts(material_name)
    material_name = check_parentheses(material_name)

    return material_name.strip()
