# coding=utf-8
import regex as re
import sympy as smp
import material_parser.core.chemical_sets as cs
import material_parser.core.constants as C


class RegExParser:
    def __init__(self):
        self._list_of_elements = cs.list_of_elements
        self._greek_symbols = "".join(C.GREEK_CHARS)
        self._comparing_symbols = "".join(C.COMPARE_SIGNS)

    def separate_phase(self, formula):
        """
        separate phase symbol part from formula
        :param formula: material string
        :return: phase symbol(s) and rest of the formula
        """
        re_phase_prefix = r"^([A-Za-z" + self._greek_symbols + r"][0-9]{0,1})\-[A-Z]\.*"
        phase = ""
        start = 0
        for m in re.finditer(re_phase_prefix, formula):
            phase = m.group(1)
            start = m.end() - 1

        return phase, formula[start:]

    def get_oxygen_deficiency(self, formula):
        """
        separate oxygen deficiency from formula
        :param formula:
        :return:
        """
        re_signs = "[" + "".join(C.SIGNS) + "]".replace("+", "\+").replace("-", "\-")
        re_oxygen_deficiency = r"O[0-9]*([±\+\-∓]{1})[a-z" + self._greek_symbols + r"]{1}[0-9]*$"

        formula_upd = formula
        oxy_def = ""
        oxy_def_sym = ""

        if len(formula) < 3:
            return formula_upd, oxy_def, oxy_def_sym

        if formula[-2:] in cs.list_of_elements_2:
            return formula_upd, oxy_def, oxy_def_sym

        for m in re.finditer(re_oxygen_deficiency, formula_upd.rstrip(")")):
            end = formula_upd[m.start():m.end()]
            splt = re.split(re_signs, end)
            oxy_def_sym = splt[-1]
            oxy_def = m.group(1)
            formula_upd = formula_upd[:m.start()] + formula_upd[m.start():].replace(end, splt[0])

        if oxy_def_sym not in C.SIGNS and oxy_def_sym == formula_upd.rstrip(")")[-1]:
            oxy_def_sym = "±"

        return formula_upd, oxy_def, oxy_def_sym

    def make_fraction_convertion(self, formula):
        r_a = r"([0-9\.]*)"
        r_b = r"(\([0-9\.]*)"
        r_x = r"([a-z]*)"
        r_s = r"([\-\+]+)"
        r_d = r"([0-9\.]*)"
        r_y = r"([a-z]+\))"
        r_c = r"(?=[/]*([0-9\.]*))"
        re_formula_fraction = r"(" + r_a + r_b + r_x + r_s + r_d + r_y + r_c + r")"
        formula_upd = formula
        for m in re.finditer(re_formula_fraction, formula_upd):
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

    def convert_weird_syntax(self, formula):
        """
        check for any weird syntax (A,B)zElxEly...
        replacing with MzElxEly... and M = [A, B]
        :param formula:
        :return:
        """
        re_weird_syntax = r"(\([A-Za-z\s]+[\/,\s]+[A-Za-z]+\))"
        variables = []
        for m in re.finditer(re_weird_syntax, formula):
            variables = re.split(r"[\/,]", m.group(0).strip('()'))
            formula = formula.replace(m.group(0), "M", 1)
        return formula, variables

    """
    formula processing: weird syntax in formulas: (A,B)zElxEly...
    """
    @property
    def re_weird_syntax(self):
        return r"(\([A-Za-z\s]+[\/,\s]+[A-Za-z]+\))"

    """
    formula processing: splitting into elements
    """
    @property
    def re_sym_dict(self):
        return r"([A-Z□]{1}[a-z]{0,1})\s*([\-\*\.\da-z" + self._greek_symbols + r"\+\/]*)"

    """
    formula processing: expression in parentheses
    """
    @property
    def re_in_parentheses(self):
        return r"\(((?>[^\(\)]+|(?R))*)\)\s*([-*\.\da-z\+/]*)"

    """
    formula processing: finsing stoichiometric variables
    """
    @property
    def re_variables(self):
        return r"[a-z" + self._greek_symbols + r"]"

    """
    additives processing
    """
    @property
    def re_additive_fraction(self):
        return  r"[\-\+:·]{0,1}\s*[0-9x\.]*\s*[vmolwt\s]*\%"

    @property
    def re_split_mixture(self):
        return r"(?<=[0-9\)])[\-⋅·∙\∗](?=[\(0-9](?!x))|" + \
               r"(?<=[A-Z])[\-⋅·∙\∗](?=[\(0-9])|" + \
               r"(?<=[A-Z\)])[\-⋅·∙\∗](?=[A-Z])|" + \
               r"(?<=[0-9\)])[\-⋅·∙\∗](?=[A-Z])" + \
               "".join([r"|(?<=" + e + r")[\-⋅·∙\∗](?=[\(0-9A-Z])" for e in self._list_of_elements]) + \
               r"|[-·]([nx0-9\.]H2O)"

    """
    mixture splitting
    """
    @property
    def re_split_mixture_prefix(self):
        return r"(^\(1\-[xyz][-xyz]*\))|(^\(100\-[xyz][\-xyz]*\))"

    @property
    def re_split_mixture_separators(self):
        return r"(.*)[\-\+·∙\∗⋅]"

    @property
    def re_split_mixture_refine(self):
        return r"(?<=[A-Z\)])[\-·∙\∗⋅](?=[xyz])|(?<=O[0-9\)]+)[\-·∙\∗⋅](?=[xyz])"

    """
    finding stoichiometric variables
    """
    @property
    def re_stoichiometric_values(self):
        return r"\s*=\s*([-]{0,1}[0-9\.\,/and\s]+)[\s\)\]\,]"

    @property
    def re_stoichiometric_range_lhs(self):
        return r"([0-9\.\s]*)\s*[<≤⩽]{0,1}\s*"

    @property
    def re_stoichiometric_range_rhs(self):
        return r"\s*[<≤⩽>]{1}\s*([0-9\.\s]+)[\s\)\]\.\,]"

    @property
    def re_stoichiometric_range_hyphen(self):
        return r"\s*=\s*([0-9\.]+)\s*[-–]\s*([0-9\.\s]+)[\s\)\]\,m\%]"

    @property
    def re_stoichiometric_range_ft(self):
        return r"[a-z\s]*from\s([0-9\./]+)\sto\s([0-9\./]+)"

    """
    finding elements variables
    """
    @property
    def re_elements_values(self):
        return r"\s*[=:]{1}\s*([A-Za-z0-9\+,\s]+)"


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

"""
stoichiometric variables
"""
re_stoichiometric_values = r"\s*=\s*([-]{0,1}[0-9\.\,/and\s]+)[\s\)\]\,]"
re_stoichiometric_range_lhs = r"([0-9\.\s]*)\s*[<≤⩽]{0,1}\s*"
re_stoichiometric_range_rhs = r"\s*[<≤⩽>]{1}\s*([0-9\.\s]+)[\s\)\]\.\,]"
re_stoichiometric_range_hyphen = r"\s*=\s*([0-9\.]+)\s*[-–]\s*([0-9\.\s]+)[\s\)\]\,m\%]"
re_stoichiometric_range_ft = r"[a-z\s]*from\s([0-9\./]+)\sto\s([0-9\./]+)"

"""
elements variables
"""
re_elements_values = r"\s*[=:]{1}\s*([A-Za-z0-9\+,\s]+)"

"""
acronyms dictionary
"""
re_capitals_no_O = "[A-NP-Z]"
