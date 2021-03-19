import abc


class PostprocessingABC(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'process_data') and
                callable(subclass.process_data) or
                NotImplemented)

    @abc.abstractmethod
    def process_data(self, chemical_structure, text_sentences=[]):
        """
        postprocess the chemical structure and find additional information in the surrounding text (if provided)
        :param chemical_structure:
        :param text_sentences:
        :return:
        """
        raise NotImplementedError