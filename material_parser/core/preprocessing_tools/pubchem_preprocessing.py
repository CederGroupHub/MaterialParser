import pubchempy as pcp
import regex as re
from material_parser.core.preprocessing_tools.preprocessing_abc import PreprocessingABC


class PubchemProcessing(PreprocessingABC):

    def process_string(self, material_string, chemical_structure):
        """
        looking up material formula in PubChem database
        THIS STEP SIGNIFICANTLY SLOWS DOWN THE OVERALL CODE PERFORMANCE
        :param material_string:
        :param chemical_structure:
        :return:
        """
        material_formula = material_string
        if re.findall("[a-z]{4,}", material_string) and not chemical_structure.material_formula:
            pcp_compounds = pcp.get_compounds(material_string, 'name')
            if len(pcp_compounds) != 0:
                material_formula = pcp_compounds[0].molecular_formula

        if material_formula != material_string:
            chemical_structure.material_name = material_string
            chemical_structure.material_formula = material_formula

        return material_formula, chemical_structure

