from .fragment_ion import FragmentIon, FragmentIonException
from .utils import AMINO_MONO_MASS, MASS_DICT


class CommonFragmentIonException(FragmentIonException):
    pass


class CommonFragmentIon(FragmentIon):
    """A derived class to encapsulate the concept of
    a common fragment ion, that is a fragment ion that
    does not possess a linker.
    """

    def __init__(
        self, pep_id, ion_type, pep_rep
    ):
        self.pep_id = pep_id
        self.ion_type = ion_type
        self.pep_rep = pep_rep
        self._integrity_check()
        self.mass = self._calc_mass()

    def __repr__(self):
        return "CommonFragmentIon: %s %s - Mass: %0.2f Da" % (
            self.get_roepstorff(),
            self.get_sequence(),
            self.mass
        )

    def _integrity_check(self):
        """
        Confirms that the:
            ion type is in "abc" or "xyz",
            pep_id is of "alpha" or "beta",
            pep_rep is a string.
        """
        if self.pep_id not in ("alpha", "beta"):
            raise CommonFragmentIonException(
                "Cannot create common fragment ion. "
                "peptide ID '%s' should be alpha or beta" % self.pep_id
            )
        if self.ion_type not in "abcxyz":
            raise CommonFragmentIonException(
                "Cannot create common frgament ion. "
                "Ion Type '%s' is not in abcxzy." % self.ion_type
            )
        aa_keys = AMINO_MONO_MASS.keys()
        for aa in self.pep_rep:
            if aa not in aa_keys:
                raise CommonFragmentIonException(
                    "Cannot create common fragment ion. "
                    "Peptide representation contains unknown amino acids"
                )

    @property
    def direction(self):
        """
        Finds the direction of the ion sequence based on the ion type.
        """
        if self.ion_type in "abc":
            return "N"
        else:
            return "C"

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
        """
        mass = 0.0
        for aa in self.pep_rep:
            mass += AMINO_MONO_MASS[aa]
        frag_dict = {
            "a": 1.0*MASS_DICT["H"] - 1.0*(MASS_DICT["C"]+MASS_DICT["O"]),
            "b": 1.0*MASS_DICT["H"],
            "c": 4.0*MASS_DICT["H"] + 1.0*MASS_DICT["N"],
            "x": 1.0*MASS_DICT["H"] + 2.0*MASS_DICT['O'] + 1.0*MASS_DICT['C'],
            "y": 3.0*MASS_DICT["H"] + MASS_DICT['O'],
            "z": 1.0*MASS_DICT['H'] + 1.0*MASS_DICT['O'] - 1.0*MASS_DICT["N"]
        }
        mass += frag_dict[self.ion_type]
        return mass

    def get_mass(self):
        """
        Return the mass of the common fragment ion.
        """
        return self.mass

    def get_roepstorff(self):
        """
        Outputs the Roepstorff nomenclature of the
        common fragment ion.

        e.g. -> "Ab1"
        """
        return "%s%s%s" % (
            "A" if self.pep_id == "alpha" else "B",
            self.ion_type,
            len(self.pep_rep)
        )

    def get_sequence(self):
        """
        Outputs the string representation as a tuple
        of the fragment ion in the correct direction.
        """
        return (self.pep_rep,)

    def get_sequence_str(self):
        """
        Outputs the string representation of the fragment
        ion in the correct direction. Useful for JSON output!
        """
        return self.pep_rep

    def to_tuple(self):
        """
        Return a tupule of form
        ("Roepstorff", mass, ("sequence",))
        """
        return (
            self.get_roepstorff(),
            self.get_mass(),
            self.get_sequence()
        )

    def ion_name(self):
        """
        Returns type of ion i.e. cross-linked, immonium etc.
        """
        return "common"
