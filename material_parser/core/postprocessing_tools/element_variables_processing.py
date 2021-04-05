import regex as re
import material_parser.core.regex_parser as rp
import material_parser.core.chemical_sets as cs
from material_parser.core.postprocessing_tools.postprocessing_abc import PostprocessingABC


class ElementVariablesProcessing(PostprocessingABC):

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
                values = self.__get_elements_from_sentence(variable, text_sentences[i].strip('., '))
                i += 1
            updated_variables[variable] = values

        chemical_structure.elements_x = updated_variables

    @staticmethod
    def __get_elements_from_sentence(var, sentence):
        """
        find elements values for var in the sentence
        :param var: <str> variable name
        :param sentence: <str> sentence to look for
        :return: <list> of <str> found values
        """
        values = re.findall(var + rp.re_elements_values, sentence)
        values = [c.rstrip("0987654321+") for v in values for c in re.split(r"[,\s]", v)
                  if c.rstrip("0987654321+") in cs.list_of_elements]

        return list(set(values))