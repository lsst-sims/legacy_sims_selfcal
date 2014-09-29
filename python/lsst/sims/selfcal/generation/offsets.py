import numpy as np
import numpy.lib.recfunctions as rfn


def BaseOffset(object):
    """Base class for how to make offset classes """
    def __init__(self, **kwargs):
        pass
    def run(self, stars, visit):
        pass

def OffsetsSys(BaseOffset):
    def __init__(self, error_sys=0.003):
        """Systematic error floor for photometry"""
        self.error_sys = error_sys
        self.newkey = 'dmag_sys'
    def run(self, stars,visits):
        nstars = np.size(stars)
        dmag = np.random.rand(nstars)*self.error_sys
        dmag = np.array(dmag, dtype=zip([self.newkey],[float]))
        stars = rfn.merge_arrays( [ stars, dmag], flatten=True, usemask=False)
        return stars
    

def OffsetSNR(BaseOffset):
    """ """
    def __init__(self, lsstFilter='r'):
        """
        Generate offsets based on the 5-sigma limiting depth of an observation and the brightness of the star.
        Note that this takes into account previous offsets that have been applied
        (so run this after things like vingetting).
        error_sys:  Systematic error floor for photometry (in mags)
        """
        self.lsstFilter=lsstFilter
        self.error_sys = error_sys
        self.newkey = 'dmag_snr'
        
    def calcMagErrors(magnitudes, lsstFilter='r', m5=24.7):
        """Right now this is assuming airmass of 1.2 and median sky brightness for all observations.  """
        xval=magnitudes * 0.0
        error_rand = magnitudes * 0.0
        magnitude_errors = magnitudes * 0.0
        error_data={} #data from Ivezic 2008 Table 2
        error_data['filters']=np.array(['u','g','r','i','z','y'])
        error_data['gamma']=np.array([0.037,0.038,0.039,0.039,0.040,0.040])
        matched=np.where(error_data['filters'] == self.lsstFilter)
        xval = np.power(10., 0.4*(magnitudes - m5))
        error_rand = np.sqrt( (0.04-error_data['gamma'][matched])*xval+error_data['gamma'][matched]*xval*xval)
        magnitude_errors = np.sqrt(self.error_sys*self.error_sys+error_rand*error_rand)
        dmag = np.random.rand(len(magnitudes))*magnitude 
        return magnitude_errors
        
    def run(self, stars, visit):
        temp_mag = stars[self.lsstFilter+'mag'].copy()
        keys = stars.dtype.names()
        dmagKeys = [key for key in keys if key[0:4]=='dmag']
        # calc what magnitude the star has when it hits the silicon.
        for key in dmagKeys:
            temp_mag += stars[key]
        
        dmag = self.calcMagErrors(temp_mag, m5 = visit['fiveSigmaDepth'] )
        dmag = np.array(dmag, dtype=zip([self.newkey],[float]))
        stars = rfn.merge_arrays( [ stars, dmag], flatten=True, usemask=False)
        return stars
    
