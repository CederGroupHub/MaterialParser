import sympy as smp
import regex as re

from sympy.abc import _clash
import constants as C
import regex_parser as rp


def simplify(value):
    """
    simplifying stoichiometric expression
    :param value: string
    :return: string
    """

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
    for p in rp.re_parentheses_open:
        new_formula = new_formula.replace(p, "(")
    for p in rp.re_parentheses_close:
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