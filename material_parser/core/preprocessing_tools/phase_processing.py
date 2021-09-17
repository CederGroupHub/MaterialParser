from material_parser.core.preprocessing_tools.preprocessing_abc import PreprocessingABC


class PhaseProcessing(PreprocessingABC):

    def __init__(self, regex_parser):
        super(PhaseProcessing, self).__init__(regex_parser)

    def process_string(self, material_string, chemical_structure):

        if chemical_structure.material_formula:
            material_string = chemical_structure.material_formula

        phase, material_string = self._re.separate_phase(material_string)

        if phase:
            chemical_structure.phase = phase

        return material_string, chemical_structure
