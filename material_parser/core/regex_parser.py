# coding=utf-8
import regex as re
import material_parser.core.chemical_sets as cs
import material_parser.core.constants as C

"""
oxygen deficiency
"""
re_signs = r"[±\+\-∓]"
re_oxygen_deficiency = r"O[0-9]*([±\+\-∓]{1})[a-z" + "".join(C.GREEK_CHARS) + r"]{1}[0-9]*$"

"""
formula fractions
"""
r_a = r"([0-9\.]*)"
r_b = r"(\([0-9\.]*)"
r_x = r"([a-z]*)"
r_s = r"([\-\+]+)"
r_d = r"([0-9\.]*)"
r_y = r"([a-z]+\))"
r_c = r"(?=[/]*([0-9\.]*))"
re_formula_fraction = r"(" + r_a + r_b + r_x + r_s + r_d + r_y + r_c + r")"

"""
weird syntax in formulas: (A,B)zElxEly...
"""
re_weird_syntax = r"(\([A-Za-z\s]+[\/,\s]+[A-Za-z]+\))"

"""
splitting into elements
"""
re_sym_dict = r"([A-Z□]{1}[a-z]{0,1})\s*([\-\*\.\da-z" + "".join(C.GREEK_CHARS) + r"\+\/]*)"

"""
expression in parentheses
"""
re_in_parentheses = r"\(((?>[^\(\)]+|(?R))*)\)\s*([-*\.\da-z\+/]*)"

"""
stoichiometric variables
"""
re_variables = r"[a-z" + "".join(C.GREEK_CHARS) + r"]"

"""
phase prefix
"""
re_phase_prefix = r"^([A-Za-z" + "".join(C.GREEK_CHARS) + r"][0-9]{0,1})\-[A-Z]\.*"

"""
fraction additives
"""
re_additive_fraction = r"[\-\+:·]{0,1}\s*[0-9x\.]*\s*[vmolwt\s]*\%"

"""
mixture splitting
"""
re_split_mixture = r"(?<=[0-9\)])[\-⋅·∙\∗](?=[\(0-9](?!x))|" + \
                   r"(?<=[A-Z])[\-⋅·∙\∗](?=[\(0-9])|" + \
                   r"(?<=[A-Z\)])[\-⋅·∙\∗](?=[A-Z])|" + \
                   r"(?<=[0-9\)])[\-⋅·∙\∗](?=[A-Z])" + \
                   "".join([r"|(?<=" + e + r")[\-⋅·∙\∗](?=[\(0-9A-Z])" for e in cs.list_of_elements]) + \
                   r"|[-·]([nx0-9\.]H2O)"
re_split_prefix = r"(^\(1\-[xyz][-xyz]*\))|(^\(100\-[xyz][\-xyz]*\))"
re_separators = r"(.*)[\-\+·∙\∗⋅]"
re_split_mixture_2 = r"(?<=[A-Z\)])[\-·∙\∗⋅](?=[xyz])|(?<=O[0-9\)]+)[\-·∙\∗⋅](?=[xyz])"

"""
string cleanup
"""
re_compare_signs = r"\({0,1}[0-9\.]*\s*[" + \
                   "".join(C.COMPARE_SIGNS) + \
                   r"]{0,1}\s*[x,y]{0,1}\s*[" + \
                   "".join(C.COMPARE_SIGNS) + \
                   r"=]\s*[0-9\.-]*\){0,1}"

re_oxygen_deficiency_replace = r"O[0-9\.]*\s*[±\+\-∓]\s*([" + r"".join(C.GREEK_CHARS) + r"]{1})"

re_tail_trash = [r"\(⩾99", r"\(99", r"\(98", r"\(90", r"\(95", r"\(96", r"\(Alfa", r"\(Aldrich", r"\(A.R.",
                 r"\(Aladdin", r"\(Sigma", r"\(A.G", r"\(Fuchen", r"\(Furuuchi", r"\(AR\)", "（x", r"\(x",
                 r"\(Acr[a-z]*", r"\(Koj", r"\(Sho", r"\(＞99"]

re_trash_syms = ["#", "$", "!", "@", "©", "®", chr(8201), "Ⓡ", "\u200b"]

re_trash_terms = ["powder", "ceramic", "rear", "earth", "micro", "nano", "coat", "crystal", "particl", "glass"]

re_parentheses_open = ["{", "["]
re_parentheses_close = ["}", "]"]
