import material_parser.core.chemical_sets as cs

from material_parser.core.chemical_structure import ChemicalStructure, Compound
from material_parser.core.default_processing import DefaultProcessing

from pprint import pprint


class MaterialParser:
    def __init__(self, options, regex_parser, preprocessings, postprocessings):
        #print ("initializing MP")
        self._options = options
        self._regex_parser = regex_parser
        self._default_processing = DefaultProcessing(regex_parser)
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
            material_string, output_structure = p(self._regex_parser).process_string(material_string,
                                                                                     output_structure)
            #print (p.__class__.__name__)
            # print ("String: {}, formula: {}, name: {}".format(material_string,
            #                                                   output_structure.material_formula,
            #                                                   output_structure.material_name))

        """
        default functionality: extraction of data from chemical formula
        """
        material_string, output_structure = self._default_processing.process_string(material_string,
                                                                                    output_structure)
        #pprint(output_structure.to_dict())

        """
        postprocessing steps
        """
        for p in self._postprocessings:
            output_structure = p(self._regex_parser).process_data(output_structure)

        output_structure.combine_formula()

        return output_structure


class MaterialParserBuilder():

    def __init__(self):
        #print ("initializing MP builder")
        self._materialParser = None
        self._file_reader = None # DefaultFileReader()
        self._regex_parser = None #DefaultRegexpParser()
        self._preprocessings = []
        self._postprocessings = []

    def add_preprocessing(self, preprocessing): # -> Builder
        self._preprocessings.append(preprocessing)
        #print ("adding another preprocessings")
        return self

    def add_postprocessing(self, postprocessing): # -> Builder
        self._postprocessings.append(postprocessing)
        #print ("adding another postprocessings")
        return self

    def set_file_reader(self, fileReader): # -> Builder
        self._file_reader = file_reader
        return self

    def set_regex_parser(self, regex_parser): # -> Builder
        self._regex_parser = regex_parser()
        return self

    def build(self, options=None): # -> MaterialParser
        #data = fileReader.read()
        #print ("building parser")
        return MaterialParser(options,
                              self._regex_parser,
                              self._preprocessings,
                              self._postprocessings)
