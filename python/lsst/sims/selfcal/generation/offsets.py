import numpy as np
import numpy.lib.recfunctions as rfn
from lsst.sims.selfcal.clouds.Arma import ArmaSf, Clouds



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

class OffsetClouds(BaseOffset):
    def __init__(self,  sampling=256, fov=3.5):
        self.fov = fov
        self.newkey = 'dmag_cloud'
        self.SF = ArmaSf()
        self.cloud = Clouds()

    def run(self, stars, visits, **kwargs):
        # XXX-Double check extinction is close to the Opsim transparency
        extinc_mags = visits['transparency']
        if extinc_mags != 0.:
            # need to decide on how to get extinc_mags from Opsim
            # Maybe push some of these params up to be setable?
            SFtheta, SFsf = self.SF.CloudSf(500., 300., 5., extinc_mags, .55)
            # Call the Clouds
            self.cloud.makeCloudImage(SFtheta,SFsf,extinc_mags, fov=self.fov)
            # Interpolate clouds to correct position.  Nearest neighbor for speed?
            nim = self.cloud.cloudimage[0,:].size
            # calc position in cloud image of each star
            starx_interp = (np.degrees(stars['x']) + self.fov/2.)*3600./ self.cloud.pixscale
            stary_interp = (np.degrees(stars['y']) + self.fov/2.)*3600./ self.cloud.pixscale

            # Round off position and make it an int
            starx_interp = np.round(starx_interp).astype(int)
            stary_interp = np.round(stary_interp).astype(int)

            # Handle any stars that are out of the field for some reason
            starx_interp[np.where(starx_interp < 0)] = 0
            starx_interp[np.where(starx_interp > nim-1)] = nim-1
            stary_interp[np.where(stary_interp < 0)] = 0
            stary_interp[np.where(stary_interp > nim-1)] = nim-1

            dmag = self.cloud.cloudimage[starx_interp,stary_interp]
        else:
            dmag = np.zeros(stars.size)
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

