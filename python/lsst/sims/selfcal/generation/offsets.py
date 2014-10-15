import numpy as np
import numpy.lib.recfunctions as rfn


class BaseOffset(object):
    """Base class for how to make offset classes """
    def __init__(self, **kwargs):
        self.newkey = 'dmag_keyname'
        pass
    def run(self, stars, visit, **kwargs):
        pass


class NoOffset(BaseOffset):
    def __init__(self):
        """ Make no changes to the mags """
        self.newkey = 'dmag_zero'
    def run(self, stars,visits, **kwargs):
        dmag = np.zeros(stars.size, dtype=zip([self.newkey],[float]))
        return dmag

class OffsetSys(BaseOffset):
    def __init__(self, error_sys=0.003):
        """Systematic error floor for photometry"""
        self.error_sys = error_sys
        self.newkey = 'dmag_sys'
    def run(self, stars,visits, **kwargs):
        nstars = np.size(stars)
        dmag = np.random.rand(nstars)*self.error_sys
        return dmag


class OffsetSNR(BaseOffset):
    """ """
    def __init__(self, lsstFilter='r'):
        """
        Generate offsets based on the 5-sigma limiting depth of an observation and the brightness of the star.
        Note that this takes into account previous offsets that have been applied
        (so run this after things like vingetting).
        """
        self.lsstFilter=lsstFilter
        self.newkey = 'dmag_snr'
        self.gamma = {'u':0.037, 'g':0.038, 'r':0.039, 'i':0.039,'z':0.040,'y':0.040 }
        self.gamma = self.gamma[lsstFilter]

    def calcMagErrors(self, magnitudes, lsstFilter='r', m5=24.7, errOnly=False):
        """Right now this is assuming airmass of 1.2 and median sky brightness for all observations.  """
        xval=magnitudes * 0.0
        error_rand = magnitudes * 0.0
        magnitude_errors = magnitudes * 0.0
        xval = np.power(10., 0.4*(magnitudes - m5))
        magnitude_errors = np.sqrt( (0.04-self.gamma)*xval + self.gamma*xval*xval)
        if errOnly:
            dmag = magnitude_errors
        else:
            dmag = np.random.rand(len(magnitudes))*magnitude_errors
        return dmag

    def run(self, stars, visit, dmags=None):
        if dmags is None:
            dmags = {}
        temp_mag = stars[self.lsstFilter+'mag'].copy()
        # calc what magnitude the star has when it hits the silicon. Thus we compute the SNR noise
        # AFTER things like cloud extinction and vingetting.

        for key in dmags.keys():
            temp_mag = temp_mag + dmags[key]
        dmag = self.calcMagErrors(temp_mag, m5 = visit['fiveSigmaDepth'] )
        return dmag

