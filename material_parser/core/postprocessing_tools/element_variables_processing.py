import regex as re
import material_parser.core.regex_parser as rp
import material_parser.core.chemical_sets as cs
from material_parser.core.postprocessing_tools.postprocessing_abc import PostprocessingABC


class ElementVariablesProcessing(PostprocessingABC):
    
    def __init__(self, regex_parser):
        super(ElementVariablesProcessing, self).__init__(regex_parser)

    def process_data(self, chemical_structure, text_sentences):
        """
        filling variables for stoichiometry and elements
        :param chemical_structure: <ChemicalStructure> see material_parser.core.chemical_structure.py
        :param text_sentences: <list> of <str>
        :return: Dict(max_value: upper limit
                      min_value: lower limit
                      values: <list> of <float> numeric values)
        """

        updated_variables = {}
        for variable in chemical_structure.elements_x:
            values = []
            i = 0
            while not values and i < len(text_sentences):
                values = self._re.get_elements_from_sentence(variable, text_sentences[i].strip('., '))
                i += 1
            updated_variables[variable] = values

        chemical_structure.elements_x = updated_variables
