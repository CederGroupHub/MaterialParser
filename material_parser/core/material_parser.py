import material_parser.core.chemical_sets as cs
from material_parser.core.chemical_structure import ChemicalStructure, Compound
from material_parser.core.default_processing import DefaultProcessing

from pprint import pprint


class MaterialParser:
    def __init__(self, data, regexParser, preprocessings, postprocessings):
        print ("initializing MP")
        self._data = data
        self._regexParser = regexParser
        self._default_processing = DefaultProcessing()
        self._preprocessings = preprocessings
        self._postprocessings = postprocessings

    def parse(self, material_string):
        """
        :param material_string:
        :return: chemical structure (see chemical_structure.py)
        """
        output_structure = ChemicalStructure(material_string)
        if not material_string:
            return output_structure

        """
        element-like material string
        """
        if material_string in cs.list_of_elements:
            return output_structure.element_structure(material_string)
        if material_string in cs.name2element:
            return output_structure.element_structure(cs.name2element[material_string])

        """
        material string is diatomic molecule
        """
        if material_string in cs.diatomic_molecules:
            return output_structure.diatomic_molecule_structure(material_string[0])

        """
        preprocessing steps
        """
        for p in self._preprocessings:
            material_string, output_structure = p.process_string(material_string, output_structure)

        """
        default functionality: extraction of data from chemical formula
        """
        material_string, output_structure = self._default_processing.process_string(material_string, output_structure)

        """
        postprocessing steps
        """
        for p in self._postprocessings:
            output_structure = p.process_data(output_structure)

        output_structure.combine_formula()

        return output_structure


class MaterialParserBuilder():

    def __init__(self):
        print ("initializing MP builder")
        self._materialParser = None
        self._fileReader = None # DefaultFileReader()
        self._regexParser = None #DefaultRegexpParser()
        self._preprocessings = []
        self._postprocessings = []

    def addPreprocessing(self, preprocessing): # -> Builder
        self._preprocessings.append(preprocessing)
        print ("adding another preprocessings")
        return self

    def addPostprocessing(self, postprocessing): # -> Builder
        self._postprocessings.append(postprocessing)
        print ("adding another postprocessings")
        return self

    def setFileReader(self, fileReader): # -> Builder
        self._fileReader = fileReader
        return self

    def build(self): # -> MaterialParser
        #data = fileReader.read()
        print ("building parser")
        data = None
        return MaterialParser(data,
                              self._regexParser,
                              self._preprocessings,
                              self._postprocessings)
