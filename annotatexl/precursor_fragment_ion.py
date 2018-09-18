from annotatexl.fragment_ion import FragmentIon, FragmentIonException
from .utils import MASS_DICT, AMINO_MONO_MASS

class PrecursorFragmentIonException(FragmentIonException):
    pass

class PrecursorFragmentIon(FragmentIon):
    """A derived class to encapsulate the concept of
    a cross-linked precursor ion, that is the full length
    cross-link.
    """

    def __init__(
        self, frag_reps,
        frag_topology=None
    ):
        self.frag_reps = frag_reps
        self.frag_topology = frag_topology
        self.mass = self._calc_mass()

    def __repr__(self):
        return "PrecursorFragmentIon: %s %s - Mass: %0.2f Da" % (
            self.get_roepstorff(),
            self.get_sequence(),
            self.mass
        )

    def _calc_mass(self):
        """
        Calculates the mass of the peptide according to the sum of the
        monoisotopic masses of the amino acods in the sequence.
        Monoisotopic masses can be found in utils.py.
        Adds mass of DSS/BS3 crosslinker - 2*H lost during the conjugation.
        """
        mass = 0
        # Sum peptide fragment amino acid masses
        for pep_rep in self.frag_reps:
            for aa in pep_rep:
                mass += AMINO_MONO_MASS[aa]

        # Modify the mass based on to include N' C' for 2 peptides
        mass += 2*(
            1.0*MASS_DICT["O"] + 2.0*MASS_DICT["H"]
        )+ 1.0*MASS_DICT["H"]
        # Modify the mass based on to include linker
        linker_mass = 138.0680796
        mass += linker_mass
        return mass

    def get_mass(self):
        """
        Return the mass of the cross-linked peptide.
        """
        return self.mass
        
    def get_roepstorff(self):
        """
        Outputs the Roepstorff nomenclature of the
        cross-linked precursor i.e. the full length

        e.g. -> "A10-B11"
        """
        return "A%s-B%s" % (
            len(self.frag_reps[0]),
            len(self.frag_reps[1])
        )

    def get_sequence(self):
        """
        Outputs the string representation as a tuple
        of the full cross-link.
        """
        return self.frag_reps

    def get_sequence_str(self):
        """
        Outputs the string representation of the full cross-link. 
        Useful for JSON output!
        """
        return "%s-%s" % self.frag_reps

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
        return "precursor"