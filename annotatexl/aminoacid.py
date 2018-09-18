from .utils import AMINO_MONO_MASS


class AminoAcid(object):
    """
    Uses Utils.py to calculate the mass of amino acid.
    """
    def __init__(self, aa_rep):
        self.aa_rep = aa_rep

    def __str__(self):
        return "AminoAcid '%s' - %sDa" % (
            self.aa_rep, self.get_mass()
        )

    def __repr__(self):
        return self.__str__()

    def get_mass(self):
        """
        Calculate the total amino acid mass based on
        atomic masses of elemental constituents.
        """
        mass = 0
        mass += AMINO_MONO_MASS[self.aa_rep]
        return mass
