from annotatexl.annotator.observed_ion import ObservedIon
from annotatexl.common_fragment_ion import CommonFragmentIon
from annotatexl.diagnostic_fragment_ion import DiagnosticFragmentIon
from annotatexl.immonium_fragment_ion import ImmoniumFragmentIon
from annotatexl.precursor_fragment_ion import PrecursorFragmentIon
from annotatexl.xl_fragment_ion import CrosslinkFragmentIon


class Annotator(object):
    """
    Annotation class that attempts to match theoretical fragment ions
    against observed fragment ions using either a Parts Per Million (ppm)
    or Dalton (Da) based tolerance window. Depending upon the input
    keyword arguments, the correct tolerance function is correctly
    specified.
    """

    def __init__(self, ppm_or_da="ppm", tolerance=10.0):
        self.ppm_or_da = ppm_or_da
        self.tolerance = tolerance
        if self.ppm_or_da == "ppm":
            self._tol_func = self._calc_tolerance_ppm
        else:
            self._tol_func = self._calc_tolerance_da

    def _calc_tolerance_ppm(self, mass):
        """
        Calculates the match mass tolerance using ppm for each observed ion.
        """
        return mass * self.tolerance / 1e6

    def _calc_tolerance_da(self):
        """
        Calculates the match +/- mass tolerance using Da for each 
        observed ion. 
        """
        return self.tolerance / 2.0

    def _sort_list_on_mass(self, ion_list):
        """
        Takes ion list and assumes each has a get mass function.
        Sorts the list by the objects mass.
        """
        return sorted(
            ion_list, key=lambda i: i.get_mass()
        )

    def annotate(self, fragment_ion_list, observed_ion_list):
        """
        Matched each observed ion in the input peak list to the list of 
        theoretical ions generated for a cross-link. 
        Uses input match tolerance parameter specified by the user to first
        calculate the +/- mass tolerance of a match. Then generates the
        mass difference between an observed ion and a theoretical ion.
        If the difference is less than the tolerance parameter the ion is 
        matched. Returns a list of ion objects containg observed ion mass and 
        theoretical mass with the calculated difference for all matches and 
        the observed mass if not matched. 
        """
        oi_sorted = self._sort_list_on_mass(observed_ion_list)
        fi_sorted = self._sort_list_on_mass(fragment_ion_list)
        matched_list = []
        for i, oi in enumerate(oi_sorted):
            eps = self._tol_func(oi.get_mass())
            matched = False
            for j, fi in enumerate(fi_sorted):
                err = abs(oi.get_mass() - fi.get_mass())
                if err <= eps:
                    matched_list.append( (oi, fi, err) )
                    matched = True
            if not matched:
                matched_list.append( (oi, None, None) )
        return matched_list

    @staticmethod
    def matched_to_dict(matched_list):
        """
        Converts the list of tuples of (ObservedIon(...), matched ion)
        into a list of dictionaries solely containing strings and
        numbers to allow straightforward JSON serialisation.
        Calls get_roepstorff, get_sequence_str, get_mass fuctions for 
        each matched ion. Adds m/z and intensity to dictionary for all
        ions in the peak list.

        Parameters
        ----------
        matched_list: list
            The list of tuples of observed/matched pairs
        """
        matched_dict_list = []
        for obs, mat, err in matched_list:
            matched = False
            if mat is not None:
                d = {
                    "matched": True,
                    "matched_roepstorff": mat.get_roepstorff(),
                    "matched_sequence": mat.get_sequence_str(),
                    "matched_mz": mat.get_mass(),   
                    "error": err
                }

                if isinstance(mat, CrosslinkFragmentIon):
                    ion_type = "crosslink"
                elif isinstance(mat, PrecursorFragmentIon):
                    ion_type = "crosslink"
                elif isinstance(mat, CommonFragmentIon):
                    ion_type = "common"
                elif isinstance(mat, ImmoniumFragmentIon):
                    ion_type = 'immonium'
                elif isinstance(mat, DiagnosticFragmentIon):
                    ion_type = 'diagnostic'
                
            else:
                d = {
                    "matched": False,
                    "matched_ion_type": "N/A",
                    "matched_roepstorff": "N/A",
                    "matched_sequence": "N/A",
                    "matched_mz": "N/A"
                }
            d["obs_mz"] = obs.get_mass()
            d["obs_int"] = obs.intensity
            matched_dict_list.append(d)
        return matched_dict_list