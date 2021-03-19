# coding=utf-8
import regex as re
import material_parser.core.regex_parser as rp
import material_parser.core.chemical_sets as cs
from material_parser.core.preprocessing_tools.preprocessing_abc import PreprocessingABC


class AdditivesProcessing(PreprocessingABC):
    __doping_terms = {"activated", "modified", "stabilized", "doped", "added"}

    def process_string(self, material_string, chemical_structure):
        """
        resolves doped part in material string
        :param material_string:
        :param chemical_structure:
        :return:
        """
        new_material_string = material_string
        additives = []

        for t in ["codoped", "co-doped"]:
            new_material_string = new_material_string.replace(t, "doped")

        """
        separate "doped with"
        """
        for r in self.__doping_terms:
            parts = [w for w in re.split(r + " with", new_material_string) if w != ""]
            if len(parts) > 1:
                new_material_string = parts[0].strip(" -+")
                additives.append(parts[1].strip())

        """
        separate element-doped prefix
        """
        for r in self.__doping_terms:
            parts = [w for w in re.split(r"(.*)[-\s]{1}" + r + " (.*)", new_material_string) if w != ""]
            if len(parts) > 1:
                new_material_string = parts.pop()
                additives.extend(parts)

        """
        separate fractions: e.g. (K0.16Na0.84)0.5Bi4.5Ti4O15+xwt.% CeO2 -> (K0.16Na0.84)0.5Bi4.5Ti4O15 and CeO2
        """
        parts = []
        if "%" in new_material_string:
            new_material_string = new_material_string.replace(".%", "%")
            parts = re.split(rp.re_additive_fraction, new_material_string)

        if len(parts) > 1:
            new_material_string = parts[0].strip(" -+")
            additives.extend(d.strip() for d in parts[1:] if d != "")

        """
        separate element(s) before/after formula: e.g. Ba5Si8O21:0.02Eu2+,xDy3+ -> Ba5Si8O21 and Eu and Dy
        """
        for part_ in new_material_string.split(":"):
            part_ = part_.strip(" ")

            part = part_
            if any(e in part for e in cs.list_of_elements_2):
                for e in cs.list_of_elements_2:
                    part = part.replace(e, "&&")

            if all(e.strip("zyx,+0987654321. ") in cs.list_of_elements_1 | {"R", "&&"}
                   for e in re.split(r"[\s,/]", part) if e != ""):
                additives.append(part_.strip(" "))
            else:
                new_material_string = part_.strip(" ")

        additives_final = [a.strip(" ") for s in additives for a in re.split(r"[\s,\-/]|and", s) if a.strip(" ") != ""]
        chemical_structure.additives = additives_final
        new_material_string = new_material_string.rstrip(r"( ,.:;-±/+")

        return new_material_string, chemical_structure
