from .peptide import Peptide


class CrosslinkException(Exception):
    pass


class Crosslink(object):
    """
    Defines a crosslinked-peptide object represention. This
    consists of 'alpha' and 'beta' peptide strings that
    represent the crosslink, with an assoiation linker
    topology.

    Example: The following crosslink has an alpha string
    given by "PEPTID", a beta string given by "HIKE" and
    a linker topology of (4, 3):

    PEPTID
       |
     HIKE

    Parameters
    ----------
    alpha_pep_rep : str
        The alpha peptide string representation
        e.g. "PEPTID"
    beta_pep_rep : str
        The beta peptide string representation
        e.g. "HIKE"
    topology : tuple
        A two-tuple containing one-indexed integers
        that represent the linker position in the
        peptide string representations
    """

    def __init__(
        self, alpha_pep_rep, beta_pep_rep, topology
    ):
        # Peptide string representations
        self.alpha_pep_rep = alpha_pep_rep
        self.beta_pep_rep = beta_pep_rep

        # Actual Peptide objects
        self.alpha_pep = Peptide(self.alpha_pep_rep)
        self.beta_pep = Peptide(self.beta_pep_rep)

        # Crosslink topology tuples
        self.topology = topology
        self.topology_zero = (
            self.topology[0] - 1, self.topology[1] - 1
        )

    @classmethod
    def from_id(cls, crosslink_id):
        """
        Allows instantiation of class using xQuest representation of 
        crosslink. "SHCIAEVEKDAIPENLPPLTADFAEDK-DVCKNYQEAK-a20-b4"
        """
        try:
            id_spl = crosslink_id.split("-")
            assert(len(id_spl) == 4)
        except Exception:
            raise ValueError(
                "Crosslink ID is not valid format. "
                "Cannot create Crosslink."
            )

        alpha = id_spl[0]
        beta = id_spl[1]
        topo = (int(id_spl[2][1:]), int(id_spl[3][1:]))
        return cls(alpha, beta, topo)

    def get_linked_amino_acid(self, peptide_id):
        """
        Identifies the amino acid involved in the crosslink for both peptides.
        """
        if peptide_id not in ("alpha", "beta"):
            raise CrosslinkException(
                "Peptidde ID '%s' is not 'alpha' or 'beta' when trying to " 
                "find linker position." % peptide_id
            )
        if peptide_id == "alpha":
            return self.alpha_pep_rep[self.topology_zero[0]]
        else:
            return self.beta_pep_rep[self.topology_zero[1]]

    def get_unique_amino_acids(self):
        """
        Returns the unique aminod acids from the crosslink not including the
        amino acids involbed in the link.
        """
        t0, t1 = self.topology_zero[0], self.topology_zero[1]
        new_alpha = self.alpha_pep_rep[:t0] + self.alpha_pep_rep[t0 + 1:] 
        new_beta = self.beta_pep_rep[:t1] + self.beta_pep_rep[t1 + 1:]
        return list(set(new_alpha + new_beta))