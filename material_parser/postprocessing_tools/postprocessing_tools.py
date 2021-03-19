import abc


class PostrocessingTools(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'process_data') and
                callable(subclass.process_data) or
                NotImplemented)

    @abc.abstractmethod
    def process_data(self, material_data):
        """Load in the data set"""
        raise NotImplementedError