# coding=utf-8
"""
SETS
"""
NUMBERS = {chr(i) for i in range(48, 58)}
OPERATIONS = {"+", "-", "/", "*"}
VARIABLES = {"x", "y", "z", "a", "b", "c", "k", "n", "m", "d"}
LATIN_CAPITAL = {chr(i) for i in range(65, 91)}
LATIN_LOWERCASE = {chr(i) for i in range(97, 123)}
PUNCTUATIONS = {".", "(", ")"}
PUNCTUATIONS_EXT = PUNCTUATIONS | {",", "-", "[", "]", "%"}
GREEK_CHARS = {chr(i) for i in range(945, 970)}
SIGNS = {"+", "-", "±", "∓"} #plus-minus, minus-plus
DEFICIENCY_CHARS = SIGNS | {chr(948), ""} #delta, empty

DASHES = {chr(8722)} | {chr(i) for i in range(8208, 8214)}
DOTS = {chr(i) for i in [42, 215, 8226, 8270, 8727, 8729, 8901, 9679, 215, 65106, 65381, 12539, 9072]}
ARROWS = {"→", "⟶", "↑", "↓", "↔", "⇌", "⇒", "⇔", "⟹"}
COMPARE_SIGNS = {"⩽", "≤", "<"}
VACANCIES = {"□"}
HYDRATE = {chr(183)}

NAME_SYMBOLS_SET = LATIN_CAPITAL | LATIN_LOWERCASE | PUNCTUATIONS_EXT | GREEK_CHARS
FORMULA_SYMBOLS_SET = NUMBERS | LATIN_CAPITAL | LATIN_LOWERCASE | OPERATIONS | PUNCTUATIONS | GREEK_CHARS | \
                      HYDRATE | VACANCIES | {","} | DEFICIENCY_CHARS
AMOUNTS_SYMBOLS_SET = NUMBERS | VARIABLES | OPERATIONS | PUNCTUATIONS

"""
REGEX STRINGS
"""
NUMBERS_STR = "0987654321+-"
ANY_LOWERCASE = "[a-z]*"
ONE_LOWERCASE = "[a-z]+"
ANY_UPPERCASE = "[A-Z]*"
ONE_UPPERCASE = "[A-Z]+"

#chr(ord("\ue5f8"))
#chr(173)
