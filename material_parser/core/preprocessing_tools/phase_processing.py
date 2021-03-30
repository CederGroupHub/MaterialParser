import regex as re
import material_parser.core.regex_parser as rp
from material_parser.core.preprocessing_tools.preprocessing_abc import PreprocessingABC
from material_parser.core.formula_processing import separate_phase


class PhaseProcessing(PreprocessingABC):
    def process_string(self, material_string, chemical_structure):

        if chemical_structure.material_formula:
            material_string = chemical_structure.material_formula

        phase, material_string = separate_phase(material_string)

        if phase:
            chemical_structure.phase = phase

        return material_string, chemical_structure
