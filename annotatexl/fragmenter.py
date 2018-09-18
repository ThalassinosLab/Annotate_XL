import itertools
import pprint

from annotatexl.common_fragment_ion import CommonFragmentIon
from annotatexl.crosslink import Crosslink
from annotatexl.immonium_fragment_ion import ImmoniumFragmentIon
from annotatexl.diagnostic_fragment_ion import DiagnosticFragmentIon
from annotatexl.precursor_fragment_ion import PrecursorFragmentIon
from annotatexl.xl_fragment_ion import CrosslinkFragmentIon


class Fragmenter(object):
    def __init__(self):
        pass

    def _n_frag_ion_strs(self, pep_rep):
        """
        Generator to create the fragment ions from the N' side
        of the peptide. E.g. "PEPTID" creates ['P', 'PE', 'PEP',
        'PEPT', 'PEPTI', 'PEPTID']
        """
        # +1 includes full length for xl frag
        for aa in range(1, len(pep_rep)+1):
            yield pep_rep[0:aa]

    def _c_frag_ion_strs(self, pep_rep):
        """
        Reverses string to get C' ions. Uses _n_frag_ion_strs to
        generate all C' ions. E.g. "PEPTID" is first reversed to
        "DITPEP" to create ['D', 'DI', 'DIT', 'DITP', 'DITPE',
        'DITPEP']
        """
        rev = "".join(list(reversed(pep_rep)))
        return self._n_frag_ion_strs(rev)

    def _n_frag_ion_strs_up_to_linker(self, pep_rep, position):
        """
        Generator to create the fragment ion strings for the
        N' common fragment ions i.e. the abc ions up to but not
        including the linker. E.g. "PEPTID, (4)" creates
        ['P'. 'PE', 'PEP'].
        """
        for i, frag in enumerate(self._n_frag_ion_strs(pep_rep)):
            if i < position:
                yield frag

    def _c_frag_ion_strs_up_to_linker(self, pep_rep, position):
        """
        Generator to create the fragment ion strings for the
        C' common fragment ions i.e. the xyz ions up to but not
        including the linker. E.g. "PEPTID" becomes "DITPEP" with
        forward position 4 (3 zero-indexed), which becomes
        zero-indexed position 2, to create ['D', 'DI'].
        """
        new_pos = len(pep_rep) - position - 1
        rev = "".join(list(reversed(pep_rep)))
        return self._n_frag_ion_strs_up_to_linker(rev, new_pos)

    def _common_frag_ion_strs_dict(self, crosslink):
        """
        Creates a dictionary of all common fragment ions up to but
        not including the linker position for the given crosslink
        and crosslink topology.
        """
        alpha = crosslink.alpha_pep_rep
        beta = crosslink.beta_pep_rep
        topo = crosslink.topology_zero
        return {
            'An': self._n_frag_ion_strs_up_to_linker(alpha, topo[0]),
            'Bn': self._n_frag_ion_strs_up_to_linker(beta, topo[1]),
            'Ac': self._c_frag_ion_strs_up_to_linker(alpha, topo[0]),
            'Bc': self._c_frag_ion_strs_up_to_linker(beta, topo[1])
        }

    def _immonium_frag_ions(self, crosslink):
        """
        Creates Immonium Ions from unique amino acids in crosslink.
        """
        unique_amino_acids = crosslink.get_unique_amino_acids()
        for amino_acid in unique_amino_acids:
            yield ImmoniumFragmentIon(amino_acid)

    def _diagnostic_frag_ions(self, crosslink):
        """
        Creates the diagnostic ions for BS3/DSS crosslinker
        """
        ions=['DI_1', 'DI_2']
        for di in ions:
            yield DiagnosticFragmentIon(di)

    def _precursor_fragment_ion(self, crosslink):
        """
        Generates ion for full cross-linked precursor
        """
        alpha = crosslink.alpha_pep_rep
        beta = crosslink.beta_pep_rep
        frag_rep = (alpha, beta)
        yield PrecursorFragmentIon(frag_rep)

    def _common_frag_ions(self, crosslink):
        """
        Returns a generator that produces every potential
        common fragment ions for both alpha and beta peptides
        from both N' and C' terminal directions, including all
        6 variations of 'abc' and 'xyz' ion types.
        """
        common_strs_dict = self._common_frag_ion_strs_dict(crosslink)
        # Create the An common fragment ions
        for frag in common_strs_dict["An"]:
            for ion_type in 'abc':
                yield CommonFragmentIon(
                    "alpha", ion_type, frag
                )
        # Create the Bn common fragment ions
        for frag in common_strs_dict["Bn"]:
            for ion_type in 'abc':
                yield CommonFragmentIon(
                    "beta", ion_type, frag
                )
        # Create the Ac common fragment ions
        for frag in common_strs_dict["Ac"]:
            for ion_type in 'xyz':
                yield CommonFragmentIon(
                    "alpha", ion_type, frag
                )
        # Create the Bc common fragment ions
        for frag in common_strs_dict["Bc"]:
            for ion_type in 'xyz':
                yield CommonFragmentIon(
                    "beta", ion_type, frag
                )

    def _n_frag_ion_strs_linked(self, pep_rep, position, complete=True):
        """
        Creates the linked N' fragment ions by taking the set of all possible
        fragments and fragment ions that don't include the linker,
        i.e. those before the amino acid at the linker position given
        by the crosslink topology.
        Removes the full precursor sequence from the final list.
        """
        linked = set(self._n_frag_ion_strs(pep_rep)) - \
            set(self._n_frag_ion_strs_up_to_linker(pep_rep, position))
        pep_list = sorted(linked, key=lambda x: len(x))
        if not complete:
            pep_list = pep_list[:-1]
        for pep in pep_list:
            yield pep

    def _c_frag_ion_strs_linked(self, pep_rep, position, complete=True):
        """
        Calculates the new position of the crosslinker for a reversed
        peptide sequence.
        Creates the linked C' fragment ions by using the _n_frag_ion_strs
        method and providing it with a set of reversed peptide sequences.
        """
        new_pos = len(pep_rep) - position - 1
        rev = "".join(list(reversed(pep_rep)))
        return self._n_frag_ion_strs_linked(rev, new_pos, complete)

    def _crosslinked_single_frag_ions_strs_dict(self, crosslink):
        """
        Cartesian product of N' and C' fragment ions for alpha and beta
        peptides. Results passed to a dictionary defining direction of
        fragment ion and whether A or B is complete.
        """
        alpha = crosslink.alpha_pep_rep
        beta = crosslink.beta_pep_rep
        topo = crosslink.topology_zero

        An_B_cp = itertools.product(
            self._n_frag_ion_strs_linked(alpha, topo[0], complete=True),
            [beta]
        )
        A_Bn_cp = itertools.product(
            [alpha],
            self._n_frag_ion_strs_linked(beta, topo[1], complete=True),
        )
        Ac_B_cp = itertools.product(
            self._c_frag_ion_strs_linked(alpha, topo[0], complete=True),
            [beta]
        )
        A_Bc_cp = itertools.product(
            [alpha],
            self._c_frag_ion_strs_linked(beta, topo[1], complete=True),
        )
        return {
            'A-Bn': [i for i in A_Bn_cp],
            'A-Bc': [i for i in A_Bc_cp],
            'An-B': [i for i in An_B_cp],
            'Ac-B': [i for i in Ac_B_cp]
        }

    def _crosslinked_single_frag_ions(self, crosslink):
        """
        Creates all ion types for the single fragmentation event crosslinks.
        """
        xl_single_frag_strs_dict = \
            self._crosslinked_single_frag_ions_strs_dict(
                crosslink
            )
        
        # Remove last item of cartesian product as it represents 
        # the precursor ion. Removed from list as list is sorted
        # Dict is not
        xl_single_frag_strs_dict['A-Bn'].pop()
        xl_single_frag_strs_dict['A-Bc'].pop()
        xl_single_frag_strs_dict['An-B'].pop()
        xl_single_frag_strs_dict['Ac-B'].pop()

        # Create the LAn crosslinked fragment ions
        for alpha_frag, beta_frag in xl_single_frag_strs_dict["An-B"]:
            for ion_type in 'abc':
                yield CrosslinkFragmentIon(
                    (alpha_frag, beta_frag),
                    (ion_type, None)
                )
        # Create the LBn crosslinked fragment ions
        for alpha_frag, beta_frag in xl_single_frag_strs_dict["A-Bn"]:
            for ion_type in 'abc':
                yield CrosslinkFragmentIon(
                    (alpha_frag, beta_frag),
                    (None, ion_type)
                )
        # Create the LAc crosslinked fragment ions
        for alpha_frag, beta_frag in xl_single_frag_strs_dict["Ac-B"]:
            for ion_type in 'xyz':
                yield CrosslinkFragmentIon(
                    (alpha_frag, beta_frag),
                    (ion_type, None)
                )
        # Create the LBc crosslinked fragment ions
        for alpha_frag, beta_frag in xl_single_frag_strs_dict["A-Bc"]:
            for ion_type in 'xyz':
                yield CrosslinkFragmentIon(
                    (alpha_frag, beta_frag),
                    (None, ion_type)
                )

    def _crosslinked_double_frag_ions_strs_dict(self, crosslink):
        """
        Cartesian product of N' and C' fragment ions for alpha and beta
        peptides. Results passed to a dictionary defining direction of
        fragment ion and whether A or B is complete.
        """
        alpha = crosslink.alpha_pep_rep
        beta = crosslink.beta_pep_rep
        topo = crosslink.topology_zero

        An_Bn_cp = itertools.product(
            self._n_frag_ion_strs_linked(alpha, topo[0], complete=False),
            self._n_frag_ion_strs_linked(beta, topo[1], complete=False)
        )
        Ac_Bn_cp = itertools.product(
            self._c_frag_ion_strs_linked(alpha, topo[0], complete=False),
            self._n_frag_ion_strs_linked(beta, topo[1], complete=False),
        )
        An_Bc_cp = itertools.product(
            self._n_frag_ion_strs_linked(alpha, topo[0], complete=False),
            self._c_frag_ion_strs_linked(beta, topo[1], complete=False)
        )
        Ac_Bc_cp = itertools.product(
            self._c_frag_ion_strs_linked(alpha, topo[0], complete=False),
            self._c_frag_ion_strs_linked(beta, topo[1], complete=False),
        )
        return {
            'An-Bn': An_Bn_cp,
            'Ac-Bn': Ac_Bn_cp,
            'An-Bc': An_Bc_cp,
            'Ac-Bc': Ac_Bc_cp
        }

    def _crosslinked_double_frag_ions(self, crosslink):
        """
        Creates all ion types for the double fragmentation event crosslinks.
        """
        xl_double_frag_strs_dict = \
            self._crosslinked_double_frag_ions_strs_dict(
                crosslink
            )
        # Create the An-Bn crosslinked fragment ions
        for frag_pair in xl_double_frag_strs_dict["An-Bn"]:
            for ion_pair in zip('abc', 'abc'):
                yield CrosslinkFragmentIon(
                    frag_pair, ion_pair
                )
        # Create the Ac-Bn crosslinked fragment ions
        for frag_pair in xl_double_frag_strs_dict["Ac-Bn"]:
            for ion_pair in zip('xyz', 'abc'):
                yield CrosslinkFragmentIon(
                    frag_pair, ion_pair
                )
        # Create the An-Bc crosslinked fragment ions
        for frag_pair in xl_double_frag_strs_dict["An-Bc"]:
            for ion_pair in zip('abc', 'xyz'):
                yield CrosslinkFragmentIon(
                    frag_pair, ion_pair
                )
        # Create the Ac-Bc crosslinked fragment ions
        for frag_pair in xl_double_frag_strs_dict["Ac-Bc"]:
            for ion_pair in zip('xyz', 'xyz'):
                yield CrosslinkFragmentIon(
                    frag_pair, ion_pair
                )

    def _fragment(self, crosslink):
        """
        Fragments the crosslink peptides to generate fragment ions series.
        """
        for frag in self._precursor_fragment_ion(crosslink):
            yield frag
        for frag in self._diagnostic_frag_ions(crosslink):
            yield frag
        for frag in self._immonium_frag_ions(crosslink):
            yield frag
        for frag in self._common_frag_ions(crosslink):
            yield frag
        for frag in self._crosslinked_single_frag_ions(crosslink):
            yield frag
        for frag in self._crosslinked_double_frag_ions(crosslink):
            yield frag

    def cid(self, crosslink):
        """
        Uses the fragment method on a crosslinked object to generate
        xl_fragment_ion and common_fragment_ion objects.
        """
        original_frags = self._fragment(crosslink)
        for frag in original_frags:
            yield frag
