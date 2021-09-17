import sympy as smp
import regex as re

from sympy.abc import _clash
import material_parser.core.constants as C
import material_parser.core.regex_parser as rp
import material_parser.core.chemical_sets as cs

re_parentheses_open = ["{", "["]
re_parentheses_close = ["}", "]"]

def simplify(value):
    """
    simplifying stoichiometric expression
    :param value: string
    :return: string
    """

    if not value:
        return value

    for l in C.GREEK_CHARS:
        _clash[l] = smp.Symbol(l)

    new_value = value
    for i, m in enumerate(re.finditer(r"(?<=[0-9])([a-z" + "".join(C.GREEK_CHARS) + "])", new_value)):
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


def simplify_str(value):
    value = smp.simplify(value)
    if value.is_Number:
        value = round(float(value), 3)
    value = str(value)

    return value


def cast_stoichiometry(value):
    if not value:
        return value
    value = float(value)
    if value == 1.0:
        return ""
    if value * 1000 % 1000 == 0.0:
        return str(int(value))

    return str(value)


def lcm(x, y):
    """This function takes two
    integers and returns the L.C.M."""

    lcm = None
    greater = x if x > y else y

    found = False
    while not found:
        if (greater % x == 0) and (greater % y == 0):
            lcm = greater
            found = True
        greater += 1

    return lcm

def parentheses_balanced(formula):
    opening_par = []

    for i, c in enumerate(formula):
        if c == '(':
            opening_par.append(i)
        if c == ')':
            if not opening_par:
                return False
            else:
                opening_par = opening_par[:-1]

    return opening_par == []


def check_parentheses(formula):
    """
    :param formula:
    :return:
    """

    if formula == "":
        return ""

    new_formula = formula
    for p in re_parentheses_open:
        new_formula = new_formula.replace(p, "(")
    for p in re_parentheses_close:
        new_formula = new_formula.replace(p, ")")

    if new_formula[0] == '(' and new_formula[-1] == ')' and parentheses_balanced(new_formula[1:-1]):
        new_formula = new_formula[1:-1]

    if new_formula[0] == '(' and parentheses_balanced(new_formula[1:]):
        new_formula = new_formula[1:]

    if new_formula[-1] == ')' and parentheses_balanced(new_formula[:-1]):
        new_formula = new_formula[:-1]

    if parentheses_balanced(new_formula):
        return new_formula

    par_open = []
    par_close = []
    for i, c in enumerate(new_formula):
        if c == '(':
            par_open.append(i)
        if c == ')':
            if not par_open:
                par_close.append(i)
            else:
                par_open = par_open[:-1]

    new_formula = '('*len(par_close) + new_formula + ")"*len(par_open)

    return new_formula


def is_int(num):
    try:
        return round(float(num), 3) == round(float(num), 0)
    except:
        return False


def is_materials_list(material_string):
    """
    indicates if the string is in form: element A, element B, ... + plural anion,
    e.g.: magnesium and lanthanum metals; magnesium, nickel and niobium (V) oxides
    :param material_string:
    :return:
    """
    if (any(a + "s" in material_string.lower() for a in cs.anions) or "metal" in material_string) and \
            any(w in material_string for w in ["and ", ",", " of "]):
        return True

    return False


def split_materials_list(material_string):
    """
    split material string into list of compounds when it's given in form several anions + cation,
    e.g.: magnesium and lanthanum metals; magnesium, nickel and niobium (V) oxides
    :param material_string: <str>
    :return: <list> of <str> - list of split chemical names
    """

    parts = [p for p in re.split(r"[\s\,]", material_string) if p != ""]

    """
    find tokens corresponding to anions, cations and valencies
    """
    anion = [(i, p[:-1]) for i, p in enumerate(parts) if p[:-1].lower() in cs.anions or p[:-1].lower() == "metal"]
    cation = [(i, p) for i, p in enumerate(parts) if p.lower() in cs.cations or p in cs.list_of_elements]
    valencies = [(i-1, p.strip("()")) for i, p in enumerate(parts) if p.strip("()") in cs.rome2num and i != 0]

    """
    only consider the situation when one common anion is given
    """
    if len(anion) != 1:
        return []

    result = []

    for c_i, c in cation:
        """
        combine each cation with anion (+ valency if given)
        """
        name = [cs.element2name[c]] if c in cs.element2name else [c.lower()]
        valency = "".join([v for v_i, v in valencies if v_i == c_i])
        if valency != "":
            name.append("(" + valency + ")")
        name.append(anion[0][1])

        """
        checking if hydrates
        """
        hydr_i = material_string.find("hydrate")
        if hydr_i > -1:
            """
            find hydrate prefix (if any)
            """
            pref = []
            while material_string[hydr_i - 1] != " " and hydr_i > 0:
                pref.append(material_string[hydr_i - 1])
                hydr_i -= 1
            pref = "".join([p for p in reversed(pref)])

            if pref not in cs.neg_prefixes:
                name.append(pref + "hydrate")
        result.append(" ".join(name))

    return result
