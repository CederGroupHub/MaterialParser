# coding=utf-8
import regex as re
import sympy as smp
import material_parser.core.chemical_sets as cs
import material_parser.core.constants as C


class RegExParser:
    def __init__(self):
        self._list_of_elements = cs.list_of_elements
        self._list_of_elements_1 = cs.list_of_elements_1
        self._list_of_elements_2 = cs.list_of_elements_2
        self._greek_symbols = "".join(C.GREEK_CHARS)
        self._comparing_symbols = "".join(C.COMPARE_SIGNS)
        self._doping_terms = {"activated", "modified", "stabilized", "doped", "added"}

    """
    PHASE PROCESSING
    """
    def separate_phase(self, formula):
        """
        separate phase symbol part from formula
        :param formula: material string
        :return: phase symbol(s) and rest of the formula
        """
        RE_PHASE_PREFIX = r"^([A-Za-z" + self._greek_symbols + r"][0-9]{0,1})\-[A-Z]\.*"
        phase = ""
        start = 0
        for m in re.finditer(RE_PHASE_PREFIX, formula):
            phase = m.group(1)
            start = m.end() - 1

        return phase, formula[start:]

    """
    FORMULA PROCESSING
    """
    def separate_oxygen_deficiency(self, formula):
        """
        separate oxygen deficiency from formula
        :param formula:
        :return:
        """
        RE_SIGNS = "[" + "".join(C.DEFICIENCY_CHARS) + "]".replace("+", "\+").replace("-", "\-")
        RE_OXY_DEF = r"O[0-9]*([±\+\-∓]{1})[a-z" + self._greek_symbols + r"]{1}[0-9]*$"

        formula_upd = formula
        oxy_def = ""
        oxy_def_sym = ""

        if len(formula) < 3:
            return formula_upd, oxy_def, oxy_def_sym

        if formula[-2:] in cs.list_of_elements_2:
            return formula_upd, oxy_def, oxy_def_sym

        for m in re.finditer(RE_OXY_DEF, formula_upd.rstrip(")")):
            end = formula_upd[m.start():m.end()]
            splt = re.split(RE_SIGNS, end)
            oxy_def_sym = splt[-1]
            oxy_def = m.group(1)
            formula_upd = formula_upd[:m.start()] + formula_upd[m.start():].replace(end, splt[0])

        if oxy_def_sym not in C.SIGNS and oxy_def_sym == formula_upd.rstrip(")")[-1]:
            oxy_def_sym = "±"

        return formula_upd, oxy_def, oxy_def_sym

    def make_fraction_convertion(self, formula):
        """
        converting fractions a(b+x)/c into (a/c*b+a/c*x) in formula
        :param formula:
        :return:
        """
        RE_A = r"([0-9\.]*)"
        RE_B = r"(\([0-9\.]*)"
        RE_X = r"([a-z]*)"
        RE_S = r"([\-\+]+)"
        RE_D = r"([0-9\.]*)"
        RE_Y = r"([a-z]+\))"
        RE_C = r"(?=[/]*([0-9\.]*))"
        RE_FORMULA_FRACTION = r"(" + RE_A + RE_B + RE_X + RE_S + RE_D + RE_Y + RE_C + r")"
        formula_upd = formula
        for m in re.finditer(RE_FORMULA_FRACTION, formula_upd):
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
        check and convert for any weird syntax (A,B)zElxEly...
        replacing with MzElxEly... and M = [A, B]
        :param formula:
        :return:
        """
        RE_WEIRD_SYNTAX = r"(\([A-Za-z\s]+[\/,\s]+[A-Za-z]+\))"
        variables = []
        for m in re.finditer(RE_WEIRD_SYNTAX, formula):
            variables = re.split(r"[\/,]", m.group(0).strip('()'))
            formula = formula.replace(m.group(0), "M", 1)
        return formula, variables

    """
    ADDITIVES PROCESSING
    """
    def separate_additives_fraction(self, formula):
        """
        separate fractions: e.g. (K0.16Na0.84)0.5Bi4.5Ti4O15+xwt.% CeO2 -> (K0.16Na0.84)0.5Bi4.5Ti4O15 and CeO2
        :param formula:
        :return:
        """
        parts = []
        additives = []
        RE_ADDITIVE_FRACTION = r"[\-\+:·]{0,1}\s*[0-9x\.]*\s*[vmolwt\s]*\%"
        if "%" in formula:
            formula = formula.replace(".%", "%")
            parts = re.split(RE_ADDITIVE_FRACTION, formula)

        if len(parts) > 1:
            formula = parts[0].strip(" -+")
            additives = [d.strip() for d in parts[1:] if d != ""]

        additives = [a.strip(" ") for s in additives for a in re.split(r"[\s,\-/]|and", s) if a.strip(" ") != ""]
        return  formula, additives

    def separate_doped_with(self, formula):
        """
        split "material doped with element(s)" into material and elements
        :param formula:
        :return:
        """
        additives = []
        for r in self._doping_terms:
            parts = [w for w in re.split(r + " with", formula) if w != ""]
            if len(parts) > 1:
                formula = parts[0].strip(" -+")
                additives.append(parts[1].strip())

        additives = [a.strip(" ") for s in additives for a in re.split(r"[\s,\-/]|and", s) if a.strip(" ") != ""]
        return formula, additives

    def separate_element_doped(self, formula):
        """
        split "element(s)-doped material" into element(s) and material
        :param formula:
        :return:
        """
        additives = []
        for r in self._doping_terms:
            parts = [w for w in re.split(r"(.*)[-\s]{1}" + r + " (.*)", formula) if w != ""]
            if len(parts) > 1:
                formula = parts.pop()
                additives.extend(parts)

        additives = [a.strip(" ") for s in additives for a in re.split(r"[\s,\-/]|and", s) if a.strip(" ") != ""]
        return formula, additives

    def separate_elements_colon_formula(self, formula):
        """
        separate element(s) before/after formula: e.g. Ba5Si8O21:0.02Eu2+,xDy3+ -> Ba5Si8O21 and Eu and Dy
        :param formula:
        :return:
        """
        additives = []
        for part_ in formula.split(":"):
            part_ = part_.strip(" ")

            part = part_
            if any(e in part for e in self._list_of_elements_2):
                for e in cs.list_of_elements_2:
                    part = part.replace(e, "&&")

            if all(e.strip("zyx,. "+C.NUMBERS_STR) in self._list_of_elements_1 | {"R", "&&"}
                   for e in re.split(r"[\s,/]", part) if e != ""):
                additives.append(part_.strip(" "))
            else:
                formula = part_.strip(" ")

        additives = [a.strip(" ") for s in additives for a in re.split(r"[\s,\-/]|and", s) if a.strip(" ") != ""]
        return formula, additives

    """
    MIXTURE PROCESSING
    """
    def split_mixture(self, formula):
        """
        split (x)compound1-(y)compound2-(z)compound2 into  [(compound1, x), (compound2, y), (compound2, z)]
        :param formula:
        :return:
        """
        RE_SPLIT_MIXTURE = r"(?<=[0-9\)])[\-⋅·∙\∗](?=[\(0-9](?!x))|" + \
                           r"(?<=[A-Z])[\-⋅·∙\∗](?=[\(0-9])|" + \
                           r"(?<=[A-Z\)])[\-⋅·∙\∗](?=[A-Z])|" + \
                           r"(?<=[0-9\)])[\-⋅·∙\∗](?=[A-Z])" + \
                           "".join([r"|(?<=" + e + r")[\-⋅·∙\∗](?=[\(0-9A-Z])" for e in self._list_of_elements]) + \
                           r"|[-·]([nx0-9\.]H2O)"
        RE_SPLIT_MIXTURE_REFINE = r"(?<=[A-Z\)])[\-·∙\∗⋅](?=[xyz])|(?<=O[0-9\)]+)[\-·∙\∗⋅](?=[xyz])"

        compounds = [p for p in re.split(RE_SPLIT_MIXTURE, formula) if p]
        if len(compounds) > 1:
            compounds = [p for part in compounds for p in re.split(RE_SPLIT_MIXTURE_REFINE, part)]

        if any(m.strip("0987654321") in self._list_of_elements for m in compounds[:-1]):
            compounds = ["".join([p + "-" for p in compounds]).rstrip("-")]

        """
        merge oxygen element if it gets split by mistake
        """
        merged_parts = [compounds[0]]
        for m in compounds[1:]:
            if re.findall("[A-Z]", m) == ["O"]:
                to_merge = merged_parts.pop() + "-" + m
                merged_parts.append(to_merge)
            else:
                merged_parts.append(m)

        return merged_parts

    def split_mixture_fractions(self, formula):
        """
        split (N-x-y)compound1+(x)compound2+(y)compound3 into [(compound1, N-x-y), (compound2, x), (compound2, y)]
        :param formula:
        :return:
        """
        RE_SPLIT_PREFIX = r"(^\(1\-[xyz][-xyz]*\))|(^\(100\-[xyz][\-xyz]*\))"
        RE_SEPARATORS = r"(.*)[\-\+·∙\∗⋅]"

        compounds = []
        pref = [s for s in re.split(RE_SPLIT_PREFIX, formula) if s]
        if len(pref) > 1:
            compound_temp = pref.pop()
            amount = pref.pop()
            variables = re.findall("[a-z]", amount)
            for v in variables:
                formula = formula.replace("(" + v + ")", v)
            compounds = []
            while variables:
                v = variables.pop()
                parts = re.findall(RE_SEPARATORS + v + "(.*)$", compound_temp)
                if parts:
                    compounds.append((parts[0][1], v))
                    compound_temp = parts[0][0]
            compounds.append((compound_temp, amount.strip("()")))

        return [c for c in reversed(compounds)]

    """
    PUBCHECM PROCESSING
    """
    def is_chemical_term(self, material_string):
        return re.findall("[a-z]{4,}", material_string) != []

    """
    ADDITIVES SUBSTITUTION
    """
    def get_additives_coefficient(self, additive):
        """
        find any stoichiometric coefficient next to the additive and split the list of additives
        e.g. 0.05Eu -> 0.05 and Eu
        :param additives: List
        :return:
        """
        r = r"^[x0-9\.]+|[x0-9\.]+$"
        coeff = re.findall(r, additive)
        element = [s for s in re.split(r, additive) if s != ""][0]
        return element, coeff

    def additive_symbolic_substitution(self, elements, coeff):
        """
        create symbolic expression of substition of additive into total composition
        :param elements: Compound.elements
        :param coeff:
        :return:
        """
        expr = "".join(["(" + v + ")+" for e, v in elements.items()]).rstrip("+")
        coeff = coeff[0] if not re.match("^[0]+[1-9]", coeff[0]) else "0." + coeff[0][1:]
        expr = expr + "+(" + coeff + ")"

        return expr, coeff

    """
    ELEMENTS VARIABLES PROCESSING
    """
    def get_elements_from_sentence(self, var, sentence):
        """
        find elements values for var in the sentence
        :param var: <str> variable name
        :param sentence: <str> sentence to look for
        :return: <list> of <str> found values
        """
        RE_ELEMENTS_VALUES = r"\s*[=:]{1}\s*([A-Za-z0-9\+,\s]+)"
        values = re.findall(var + RE_ELEMENTS_VALUES, sentence)
        values = [c.rstrip(C.NUMBERS_STR) for v in values for c in re.split(r"[,\s]", v)
                  if c.rstrip(C.NUMBERS_STR) in self._list_of_elements]

        return list(set(values))

    """
    formula processing: finding stoichiometric variables
    """
    @property
    def re_variables(self):
        return r"[a-z" + self._greek_symbols + r"]"

    """
    STOICHIOMETRIC VARIABLES PROCESSING
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
stoichiometric variables
"""
re_variables = r"[a-z" + "".join(C.GREEK_CHARS) + r"]"

"""
acronyms dictionary
"""
re_capitals_no_O = "[A-NP-Z]"
