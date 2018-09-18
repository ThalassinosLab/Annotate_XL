class ObservedIon(object):
    """
    Creates observed ion objects for a peak list input.
    Peak list must be a CSV with two columns titled m/z intensity.
    Returns the m/z and intensity value of each recorded ion observation.
    """

    def __init__(self, mz, intensity):
        self.mz = mz
        self.intensity = intensity

    def __str__(self):
        return "ObservedIon(%0.2f, %0.2f)" % (
            self.mz, self.intensity
        )

    def __repr__(self):
        return self.__str__()

    def get_mass(self):
        """
        Returns m/z value for each ion in the input peak list
        """
        return self.mz

    def get_intensity(self):
        """
        Returns the intensity value for each ion in the input peak list
        """
        return self.intensity