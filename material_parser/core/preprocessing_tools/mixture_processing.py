import regex as re
import material_parser.core.regex_parser as rp
import material_parser.core.chemical_sets as cs
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
            # try:
            #     f = smp.simplify(f)
            #     if f.is_Number:
            #         f = round(float(f), 3)
            #     f = str(f)
            # except:
            #     f = "1"
            f = simplify(f)
            output_compounds.append((m, f))

        if output_compounds:
            chemical_structure.material_formula = material_string

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
        pref = [s for s in re.split(rp.re_split_prefix, material_name) if s]
        if len(pref) > 1:
            material_name_temp = pref.pop()
            amount = pref.pop()
            variables = re.findall("[a-z]", amount)
            for v in variables:
                material_name = material_name.replace("(" + v + ")", v)
            compounds = []
            while variables:
                v = variables.pop()
                parts = re.findall(rp.re_separators + v + "(.*)$", material_name_temp)
                if parts:
                    compounds.append((parts[0][1], v))
                    material_name_temp = parts[0][0]
            compounds.append((material_name_temp, amount.strip("()")))
            return [c for c in reversed(compounds)]

        """
        general case xA-yB-zC
        """
        parts = [p for p in re.split(rp.re_split_mixture, material_name) if p]
        #print("-->", material_name, parts)

        parts_upd = [p for part in parts for p in re.split(rp.re_split_mixture_2, part)] if len(parts) > 1 else parts

        if any(m.strip("0987654321") in cs.list_of_elements for m in parts_upd[:-1]):
            parts_upd = ["".join([p + "-" for p in parts_upd]).rstrip("-")]

        merged_parts = [parts_upd[0]]
        for m in parts_upd[1:]:
            if re.findall("[A-Z]", m) == ["O"]:
                to_merge = merged_parts.pop() + "-" + m
                merged_parts.append(to_merge)
            else:
                merged_parts.append(m)

        composition = []
        for m in merged_parts:
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
                composition.append((m, fraction))
        return composition