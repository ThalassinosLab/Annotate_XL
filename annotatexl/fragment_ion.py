from abc import ABCMeta, abstractmethod


class FragmentIonException(Exception):
    pass


class FragmentIon(object):
    """This abstract base case specifies an interface
    to be used by all derived subclasses that are of
    type FragmentIon.

    In particular it must expose get_mass(), get_roepstorff()
    get_sequence() and to_tuple() methods. This ensures that
    any client annotation module is clear on the utilsied
    interface.
    """

    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def get_mass(self):
        raise NotImplementedError(
            "Should implement get_mass()"
        )

    @abstractmethod
    def get_roepstorff(self):
        raise NotImplementedError(
            "Should implement get_roepstorff()"
        )

    @abstractmethod
    def get_sequence(self):
        raise NotImplementedError(
            "Should implement get_sequence()"
        )

    @abstractmethod
    def to_tuple(self):
        raise NotImplementedError(
            "Should implement to_tuple()"
        )

    @abstractmethod
    def ion_name(self):
        raise NotImplementedError(
            "Should implement ion_name()"
        )