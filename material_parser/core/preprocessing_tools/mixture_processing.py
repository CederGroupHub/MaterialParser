import regex as re
from material_parser.core.preprocessing_tools.preprocessing_abc import PreprocessingABC
from material_parser.core.utils import check_parentheses, simplify
from material_parser.core.chemical_structure import Compound

from pprint import pprint


class MixtureProcessing(PreprocessingABC):

    def __init__(self, regex_parser):
        super(MixtureProcessing, self).__init__(regex_parser)

    def process_string(self, material_string, chemical_structure):
        """
        splitting the complex mixture/composite/alloy formula into compounds
        I did my best trying to make it to work for as many cases as possible,
        but there are way to many variations of notations
        this function is a pain and often leads to fails in compositions extraction
        :param material_string:
        :param chemical_structure:
        :return:
        """
        if not material_string:
            return material_string, chemical_structure

        if chemical_structure.material_formula:
            material_string = chemical_structure.material_formula

        split = self.__split_formula(material_string)
        l = 0
        while len(split) != l:
            l = len(split)
            split = [p for s in split for p in self.__split_formula(s[0], s[1])]

        output_compounds = []
        for m, f in split:
            f = simplify(f)
            output_compounds.append((m, f))

        if output_compounds:
            chemical_structure.material_formula = material_string

        """
        create composition object for output from compounds and fraction
        """
        composition = []
        formula = chemical_structure.material_formula
        for m, f in output_compounds:
            if "H2O" not in m \
                    and len([(m,f) for m,f in output_compounds if "H2O" not in m]) == 1 \
                    and f not in ["1", "1.0"]:
                material_string = re.sub("^" + f, "", material_string)
                formula = re.sub("^" + f, "", formula)
                m = re.sub("^" + f, "", m)
                f = "1"
            composition.append({"formula": check_parentheses(m),
                                "amount": f})

        chemical_structure.composition = composition
        chemical_structure.material_formula = formula

        return chemical_structure.material_formula, chemical_structure


    def __split_formula(self, material_name_, init_fraction="1"):
        """
        recursively splits complex mixture/composite/alloy formula into constituting compounds
        :param material_name_: str
        :param init_fraction:
        :return: list of tuples (compound, fraction/amount)
        """
        material_name = material_name_.replace(" ", "")

        """
        cases: (N-x-y)compound1+(x)compound2+(y)compound3
        """
        compounds = self._re.split_mixture_fractions(material_name)
        if compounds:
            return compounds

        """
        general case xA-yB-zC
        """
        compounds = self._re.split_mixture(material_name)
        return self.__prettify_compound_fraction(compounds, init_fraction)

    def __prettify_compound_fraction(self, compounds, init_fraction):
        compound_fraction = []
        for m in compounds:
            fraction = ""
            i = 0
            while i < len(m) and not m[i].isupper():
                fraction = fraction + m[i]
                i += 1
            fraction = fraction.strip("()")
            if fraction == "":
                fraction = "1"
            else:
                m = m[i:]

            fraction = "(" + fraction + ")*(" + init_fraction + ")"

            if m != "":
                compound_fraction.append((m, fraction))
        return compound_fraction