import abc


class PreprocessingABC(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'process_string') and
                callable(subclass.process_string) and
                hasattr(subclass, '_re') or
                NotImplemented)

    @abc.abstractmethod
    def __init__(self, regex_parser):
        self._re = regex_parser

    @abc.abstractmethod
    def process_string(self, material_str, chemical_structure):
        """
        preprocesses material string and extracts chemical data from there
        :param material_str: material string
        :param chemical_structure: output data structure (see material_parser/chemical_structure.py)
        :return: updated materials string, updated chemical structure
        """
        raise NotImplementedError