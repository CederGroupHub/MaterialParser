import abc


class PreprocessingTools(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'process_string') and
                callable(subclass.process_string) or
                NotImplemented)

    @abc.abstractmethod
    def process_string(self, material_str):
        """Load in the data set"""
        raise NotImplementedError