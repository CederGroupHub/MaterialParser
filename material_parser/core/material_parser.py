from chemical_structure import ChemicalStructure, Compound
from formula_processing import Formula2Composition
from default_processing import DefaultProcessing
import chemical_sets as cs

from pprint import pprint


class MaterialParser:
    def __init__(self, data, regexParser, preprocessings, postprocessings):
        print ("initializing MP")
        self.__data = data
        self.__regexParser = regexParser
        self.__default_processing = DefaultProcessing()
        self.__preprocessings = preprocessings
        self.__postprocessings = postprocessings

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
        for p in self.__preprocessings:
            material_string, output_structure = p.process_string(material_string, output_structure)

        """
        default functionality: extraction of data from chemical formula
        """
        material_string, output_structure = self.__default_processing.process_string(material_string, output_structure)

        return output_structure


class MaterialParserBuilder():

    def __init__(self):
        print ("initializing MP builder")
        self.__materialParser = None
        self.__fileReader = None # DefaultFileReader()
        self.__regexParser = None #DefaultRegexpParser()
        self.__preprocessings = []
        self.__postprocessings = []

    def addPreprocessing(self, preprocessing): # -> Builder
        self.__preprocessings.append(preprocessing)
        print ("adding another preprocessings")
        return self

    def addPostprocessing(self, postprocessing): # -> Builder
        self.__postprocessings.append(postprocessing)
        print ("adding another postprocessings")
        return self

    def setFileReader(self, fileReader): # -> Builder
        self.__fileReader = fileReader
        return self

    def build(self): # -> MaterialParser
        #data = fileReader.read()
        print ("building parser")
        data = None
        return MaterialParser(data, self.__regexParser, self.__preprocessings, self.__postprocessings)
