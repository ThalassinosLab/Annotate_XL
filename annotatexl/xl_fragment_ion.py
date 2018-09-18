from .fragment_ion import FragmentIon, FragmentIonException
from .utils import MASS_DICT, AMINO_MONO_MASS


class CrosslinkFragmentIonException(FragmentIonException):
    pass


class CrosslinkFragmentIon(FragmentIon):
    """A derived class to encapsulate the concept of
    a crosslinked fragment ion, that is a fragment ion that
    possess a linker.
    """

    def __init__(
        self, frag_reps, ion_types,
        frag_topology=None
    ):
        self.frag_reps = frag_reps
        self.ion_types = ion_types
        self.frag_topology = frag_topology
        self.mass = self._calc_mass()

    def __repr__(self):
        return "CrosslinkFragmentIon: %s %s - Mass: %0.2f Da" % (
            self.get_roepstorff(),
            self.get_sequence(),
            self.mass
        )

    def _calc_mass(self):
        """
        Calculates the mass of the peptide according to the sum of the
        monoisotopic masses of the amino acods in the sequence.
        Monoisotopic masses can be found in utils.py.
        For each ion type the in "abc" and "xyz" the calculated mass is
        adjusted to match the elemental loss due to differences in bond
        scission and presence of mobile proton.

        For example:
        "abc" ions have an additional H+, "a" have lost C=O, "c" have gained
        NH.
        "xyz" have an additional H+ and H2O, "x" have gained C=O, "z" have
        lost NH.

        Adds mass of DSS/BS3 crosslinker - 2*H lost during the conjugation.
        """
        mass = 0.0
        frag_dict = {
            "a": 1.0*MASS_DICT["H"] - 1.0*(MASS_DICT["C"]+MASS_DICT["O"]),
            "b": 1.0*MASS_DICT["H"],
            "c": 4.0*MASS_DICT["H"] + 1.0*MASS_DICT["N"],
            "x": 1.0*MASS_DICT["H"] + 2.0*MASS_DICT['O'] + 1.0*MASS_DICT['C'],
            "y": 3.0*MASS_DICT["H"] + MASS_DICT['O'],
            "z": 1.0*MASS_DICT['H'] + 1.0*MASS_DICT['O'] - 1.0*MASS_DICT["N"]
        }

        # Sum peptide fragment amino acid masses
        for pep_rep in self.frag_reps:
            for aa in pep_rep:
                mass += AMINO_MONO_MASS[aa]

        # Modify the mass based on ion type presence
        for ion_type in self.ion_types:
            if ion_type is not None:
                mass += frag_dict[ion_type]
            else:
                mass += 1.0*MASS_DICT["O"] + 2.0*MASS_DICT["H"]
        linker_mass = 138.0680796
        mass += linker_mass
        return mass

    def get_mass(self):
        """
        Return the mass of the crosslinked fragment ion.
        """
        return self.mass

    def get_roepstorff(self):
        """
        Outputs the Roepstorff nomenclature of the
        crosslinked fragment ion.

        e.g. -> "Ab1-By2"
        """
        return "A%s%s-B%s%s" % (
            self.ion_types[0] if self.ion_types[0] is not None else "",
            len(self.frag_reps[0]),
            self.ion_types[1] if self.ion_types[1] is not None else "",
            len(self.frag_reps[1])
        )

    def get_sequence(self):
        """
        Outputs the string representation as a tuple
        of the fragment ion in the correct direction.
        """
        return self.frag_reps

    def get_sequence_str(self):
        """
        Outputs the string representation of the fragment
        ion in the correct direction. Useful for JSON output!
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
        return "crosslink"
