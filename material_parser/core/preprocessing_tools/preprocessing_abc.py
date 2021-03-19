import abc


class PreprocessingABC(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'process_string') and
                callable(subclass.process_string) or
                NotImplemented)

    @abc.abstractmethod
    def process_string(self, material_str, chemical_structure):
        """
        preprocesses material string and extracts chemical data from there
        :param material_str: material string
        :param chemical_structure: output data structure (see material_parser/chemical_structure.py)
        :return: updated materials string, updated chemical structure
        """
        raise NotImplementedError