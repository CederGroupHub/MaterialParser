# coding=utf-8
from material_parser.core.preprocessing_tools.preprocessing_abc import PreprocessingABC


class AdditivesProcessing(PreprocessingABC):

    def __init__(self, regex_parser):
        super(AdditivesProcessing, self).__init__(regex_parser)

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
        new_material_string, parts = self._re.separate_doped_with(new_material_string)
        additives.extend(parts)

        """
        separate element-doped prefix
        """
        new_material_string, parts = self._re.separate_element_doped(new_material_string)
        additives.extend(parts)

        """
        separate fractions: e.g. (K0.16Na0.84)0.5Bi4.5Ti4O15+xwt.% CeO2 -> (K0.16Na0.84)0.5Bi4.5Ti4O15 and CeO2
        """
        new_material_string, parts = self._re.separate_additives_fraction(new_material_string)
        additives.extend(parts)

        """
        separate element(s) before/after formula: e.g. Ba5Si8O21:0.02Eu2+,xDy3+ -> Ba5Si8O21 and Eu and Dy
        """
        new_material_string, parts = self._re.separate_elements_colon_formula(new_material_string)
        additives.extend(parts)

        chemical_structure.additives = additives
        new_material_string = new_material_string.rstrip(r"( ,.:;-±/+")

        return new_material_string, chemical_structure
