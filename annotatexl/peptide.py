from annotatexl.aminoacid import AminoAcid


class Peptide(object):
    def __init__(self, pep_rep):
        self.pep_rep = pep_rep
        self.aa_list = self._create_aa_list()

    def _create_aa_list(self):
        """
        Creates a list of all amino acids in a peptide
        """
        aa_list = [
            AminoAcid(aa_rep)
            for aa_rep in self.pep_rep
        ]
        return aa_list

    def _get_backbone_mass(self):
        """
        Calculates the Peptide 'backbone' mass without
        the end terminal ion masses.
        """
        backbone_mass = sum([
            aa.get_mass() for aa in self.aa_list
        ])
        return backbone_mass

    def create_abc_ions(self):
        """
        Generate peptide fragment ions from N'
        """
        for aa in range(1, len(self.pep_rep)+1):
            yield self.pep_rep[0:aa]

    def create_xyz_ions(self):
        """
        Reverse pep to generate ions from C'
        """
        rev = ''.join(list(reversed(self.pep_rep)))
        for aa in range(1, len(rev)+1):
            yield rev[0:aa]

    def get_mass(self):
        """
        Calculate the total mass of the peptide.
        """
        return self._get_backbone_mass()