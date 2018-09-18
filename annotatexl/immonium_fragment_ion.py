from annotatexl.fragment_ion import FragmentIon, FragmentIonException
from annotatexl.utils import IMMON_MASS


class ImmoniumFragmentIonException(FragmentIonException):
    pass


class ImmoniumFragmentIon(FragmentIon):
    """
    Creates Immonium Ion enteries based on the Immonium Ion mass dictionary
    in utils.py. These ions require their own class as they do not require
    the same level of complexty as linear and cross-linked fragment ions.
    """

    def __init__(self, amino_rep):
        self.amino_rep = amino_rep
        self.mass = self._calc_mass()

    def __repr__(self):
        return "ImmoniumFragmentIon: %s %s - Mass: %0.2f Da" % (
            self.get_roepstorff(),
            self.get_sequence(),
            self.mass
        )

    def _calc_mass(self):
        try:
            return IMMON_MASS[self.amino_rep]
        except KeyError:
            raise ImmoniumFragmentIonException(
                "Amino Acid representation '{}' does not exist in "
                "Immonium Mass Dictionary when attempting to access "
                "get_mass method for immonium ion.".format(self.amino_rep)
            )

    def get_mass(self):
        """
        Returns immonium ion mass from Utils.py
        """
        return self.mass
        
    def get_roepstorff(self):
        """
        Defines nomenclature for immonium ion to be IM_n where n is the 
        single letter amino acid code.
        """
        return "IM_%s" % self.amino_rep

    def get_sequence(self):
        """
        Returns single letter amino acid code
        """
        return self.amino_rep

    def get_sequence_str(self):
        """
        Outputs the string representation of the immonium ion. 
        Useful for JSON output!
        """
        return self.get_sequence()

    def to_tuple(self):
        """
        Returns a tupule of form
        ("Roepstorff", mass, ("sequence1","sequence2"))
        """
        return (
            self.get_roepstorff(),
            self.get_mass(),
            self.get_sequence()
        )

    def ion_name(self):
        """
        Returns ion type i.e. cross-linked, immonium etc.
        """
        return "immonium"