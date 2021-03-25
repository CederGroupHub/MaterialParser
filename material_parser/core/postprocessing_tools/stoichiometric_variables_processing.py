import regex as re
import material_parser.core.regex_parser as rp
from material_parser.core.postprocessing_tools.postprocessing_abc import PostprocessingABC
from material_parser.core.utils import simplify


class StoichiometricVariablesProcessing(PostprocessingABC):

    def process_data(self, chemical_structure, text_sentences):
        """
        filling variables for stoichiometry
        :param chemical_structure: <ChemicalStructure> see material_parser.core.chemical_structure.py
        :param text_sentences: <list> of <str>
        :return: Dict(max_value: upper limit
                      min_value: lower limit
                      values: <list> of <float> numeric values)
        """

        updated_variables = {}
        for variable in chemical_structure.amounts_x:
            values = {"values": [],
                      "max_value": None,
                      "min_value": None}
            i = 0
            while not (len(values['values']) or values['max_value']) and i < len(text_sentences):
                values = self.__get_values_from_sentence(var, text_sentences[i])
                i += 1
            updated_variables[variable] = values

        chemical_structure.amounts_x = updated_variables

    def __get_values_from_sentence(self, var, sentence):
        """
        find numeric values of var in sentence
        :param var: <str> variable name
        :param sentence: <str> sentence to look for
        :return: <dict>: max_value: upper limit
                        min_value: lower limit
                        values: <list> of <float> numeric values
        """
        values = {"values": [],
                  "max_value": None,
                  "min_value": None}

        """
        considering 4 cases of mentioning the values of variable in the text
        1. list of discrete values
        2. range x = ... - ...
        3. range ... < x < ...
        4. x from ... to ...
        """
        regs = [(var + rp.re_stoichiometric_values, "discrete"),
                (var + rp.re_stoichiometric_range_hyphen, "range"),
                (rp.re_stoichiometric_range_lhs + var + rp.re_stoichiometric_range_rhs, "range"),
                (var + rp.re_stoichiometric_range_ft, "range")]

        for r, m in regs:
            if values["values"] == [] and values["max_value"] is None:
                r_res = re.findall(r, sentence.replace(" - ", "-"))
                values = self.__get_values(r_res, m)

        return values

    def __get_values(self, string, mode):
        if not string:
            return {"values": [],
                    "max_value": None,
                    "min_value": None}

        # given range
        if mode == "range":
            return self.__get_range_values(string)

        # given list
        if mode == "discrete":
            return self.__get_discrete_values(string)

    @staticmethod
    def __get_range_values(string):
        min_value, max_value = string[0]
        max_value = max_value.rstrip("., ")
        min_value = min_value.rstrip("., ")
        max_value = re.sub("[a-z]*", "", max_value)
        min_value = re.sub("[a-z]*", "", min_value)
        try:
            max_value = round(float(simplify(max_value)), 4)
            min_value = round(float(simplify(min_value)), 4) if min_value else None
        except Exception as ex:
            max_value = None
            min_value = None
            template = "An exception of type {0} occurred when use simplify(). Arguments:\n{1!r}."
            message = template.format(type(ex).__name__, ex.args)
            print(message)
        values = []

        return {"values": values,
                "max_value": max_value,
                "min_value": min_value}

    @staticmethod
    def __get_discrete_values(string):
        values = re.split(r"[,\s]", re.sub("[a-z]+", "", string[0]))
        try:
            values = [round(float(simplify(c.rstrip("., "))), 4) for c in values if c.rstrip("., ") not in ["", "and"]]
            max_value = max(values) if values else None
            min_value = min(values) if len(values) > 1 else None
        except Exception as ex:
            values = []
            max_value = None
            min_value = None
            template = "An exception of type {0} occurred when use simplify(). Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message)

        return {"values": values,
                "max_value": max_value,
                "min_value": min_value}