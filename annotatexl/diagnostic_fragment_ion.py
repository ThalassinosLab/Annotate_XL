from annotatexl.fragment_ion import FragmentIon, FragmentIonException
from annotatexl.utils import DIAG_IONS

class DiagnosticFragmentIonException(FragmentIonException):
    pass

class DiagnosticFragmentIon(FragmentIon):
    """
    Creates Diagnostic Ion enteries based on the Immonium Ion mass dictionary
    in utils.py. These ions require their own class as they do not require
    the same level of complexty as common and cross-linked fragment ions.
    """

    def __init__(self, diag_ion_rep):
        self.diag_ion_rep = diag_ion_rep
        self.mass = self._calc_mass()

    def __repr__(self):
        return "DiagnosticFragmentIon: %s %s - Mass: %0.2f Da" % (
            self.get_roepstorff(),
            self.get_sequence(),
            self.mass
        )

    def _calc_mass(self):
        try:
            return DIAG_IONS[self.diag_ion_rep]
        except KeyError:
            raise DiagnosticFragmentIonException(
                "Diagnostic ion representation '{}' does not exist in "
                "Diagnostic ion Dictionary when attempting to access "
                "get_mass method for diagnostic ion.".format(self.diag_ion_rep)
            )

    def get_mass(self):
        """
        Return the mass of the diagnostic ion from utils.py.
        """
        return self.mass
        
    def get_roepstorff(self):
        """
        Outputs the diagnostic ion type from utils.py.

        e.g. -> "A5-B9"
        """
        return self.diag_ion_rep

    def get_sequence(self):
        """
        Represents cross-linker fragment so sequence is replaced with a 
        description for the ion.
        """
        return "LinkerIon"

    def get_sequence_str(self):
        """
        Outputs the string representation of the linker diagnostic
        ion in the correct direction. Useful for JSON output!
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
        return "diagnostic"