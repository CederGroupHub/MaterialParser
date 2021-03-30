# coding=utf-8
NUMBERS = {chr(i) for i in range(48, 58)}
OPERATIONS = {"+", "-", "/", "*"}
VARIABLES = {"x", "y", "z", "a", "b", "c", "k", "n", "m", "d"}
LATIN_CAPITAL = {chr(i) for i in range(65, 91)}
LATIN_LOWERCASE = {chr(i) for i in range(97, 123)}
PUNCTUATIONS = {".", "(", ")"}
PUNCTUATIONS_EXT = PUNCTUATIONS | {",", "-", "[", "]", "%"}
GREEK_CHARS = {chr(i) for i in range(945, 970)}
SIGNS = {"+", "-", chr(177), chr(8723)} #plus-minus, minus-plus
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

#chr(ord("\ue5f8"))
#chr(173)
