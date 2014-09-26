"""
modtranSequence.py

Class used to generate sequences of atmospheric parameters for input to MODTRAN.

Calling method :
  seq = AtmosphereSequences('Opsim3.61.02.sel_parmlist.dat')
  seq.generate_parameters()

  modtranDictList = seq.modtran_visits
  modtranParameters = modtranDictList[i]

To Lynne:
---------
NEED TO HAVE THE ENVIRONNEMENT VARIABLE DEFINED : LSST_POINTING_DIR

For the moment, the way the code works is
  - I give to the makeatmos (e.g. modtranAtmos.Atmosphere) the mjd boundaries
of the opsim file and the number of lines of the catalog
  - Then I compute the values for each pointing by calling the appropriate
function of the Atmosphere class

6/29/12
- Possibility to input a list of opsim dictionary visits instead of an opsim_file
- checkOpsimData() to check whether the input data has no bug
- changeOpsimData(filename or data) to change input opsim data without creating a new class

10/17/12
- ozone simulator implemented
- seed implemented for the simulations
- bug fixes
"""

import numpy as np
import os

import modtranAtmos
import modtranTools

rad2dec = 180. / np.pi
deg2rad = 1./rad2dec

_opsim_keys = ['obsHistID', 'expMJD', 'fieldRA', 'fieldDEC']
#_modtran_keys = ['O3', 'H2O', 'IHAZE', 'ISEAS', 'IVULC', 'ICLD', 'IVSA', 'VIS',
#               'ZAER11', 'ZAER12', 'SCALE1', 'ZAER21', 'ZAER22', 'SCALE2']
_default_seed = 10


class AtmosphereSequence(object):
    """This class is used to generate sequences of atmospheric parameters,
    matching long-term weather trends, and accounting for the airmass/pointing
    location of each observation. These atmospheric parameters are the
    basic inputs for MODTRAN, so can be used to generate an atmospheric
    transmission curve.
    The class generates these parameters for a series of opsim pointings.
    """
    def __init__(self, opsim_filename=None, opsim_data=None):
        """Instantiate an AtmosphereSequences object.
        """
        self.opsim_filename = opsim_filename
        self.opsim_data = opsim_data
        self.opsim_visits = []
        self.modtran_visits = []

        self._initialized_visits = False
        self._initialized_sequence = False

        """
        # so, for example ..
        self.modtran_visits[i] = {}
        for key in self.modtran_keys:
            self.modtran_visits[i][key] = value
        """

    def readOpsimData(self):
        """Read the pointing file and return the content as a list
        of dictionaries called opsim_visits. The dictionary keys are
        stored in class variable_opsim_keys.
        """
        if self.opsim_data:
            good, string = self.checkOpsimData()
            if good:
                self.opsim_visits = self.opsim_data
                self._initialized_visits = True
            else:
                raise Exception(string)
        elif self.opsim_filename:
            dataDir = os.getenv('LSST_POINTING_DIR')
            if not dataDir:
                raise Exception('LSST_POINTING_DIR env not set')
            opsimfile = os.path.join(dataDir, self.opsim_filename)
            # Read the file, store the info.
            with open(opsimfile, 'r') as opsim:
                for visit in opsim:
                    if visit.startswith('o'):
                        print visit
                        continue
                    data = visit.strip().split()
                    visitDict = {}
                    visitDict[_opsim_keys[0]] = int(data[0])
                    for i in xrange(1, len(data)):
                        visitDict[_opsim_keys[i]] = float(data[i])
                    self.opsim_visits.append(visitDict)
            self._initialized_visits = True
        else:
            raise Exception('No data specified')

    def checkOpsimData(self, otherdata=None):
        """Reads thru the data to see if there's any bug
        or missing key(s)
        """
        data2check = self.opsim_data
        if otherdata:
            data2check = otherdata
        if type(data2check) != list:
            if type(data2check) == dict:
                if set(self._opsim_keys).issubset(set(data2check.keys())):
                    return True, 'good'
                else:
                    return False, 'Data is not complete'
            else:
                return False, 'Data is not a list nor a dictionary'
        else:
            nodict = 0
            baddict = 0
            for data in data2check:
                if type(data) != dict:
                    nodict += 1
                else:
                    if not set(self._opsim_keys).issubset(set(data.keys())):
                        baddict += 1
            if (nodict == 0 and baddict == 0):
                return True, 'good'
            else:
                return (False, '%d non dictionary element(s) and %s uncomplete dictionary(ies) in data list' % (nodict, baddict))

    def changeOpsimData(self, opsim_filename=None, opsim_data=None):
        """Input a new opsim data file or data list into the sequence
        """
        self.opsim_filename = opsim_filename
        self.opsim_data = opsim_data

        self._initialized_visits = False
        self._initialized_sequence = False

    def _initPointingSeq(self):
        """Initialize input parameters for the atmospheric simulation
        """
        if not self._initialized_visits:
            self.readOpsimData()
        self.npoints = len(self.opsim_visits)
        self.mjds = self.opsim_visits[0]['expMJD']
        self.mjde = self.opsim_visits[-1]['expMJD']

        self._initialized_sequence = True

    def generate_parameters(self, seed=_default_seed):
        """Generate the atmospheric parameters over time.
        Returns a list of dictionaries containing the modtran
        information for each opsim visit (in the same order).
        """
        if not self._initialized_sequence:
            self._initPointingSeq()
        self.atmos = modtranAtmos.Atmosphere(
            self.mjds, self.mjde, self.npoints, seed)
        self.atmos.init_main_parameters()
        for opsim_dict, idx in zip(self.opsim_visits,
                                   xrange(len(self.opsim_visits))):
            modtran_dict = self.fillModtranDictionary(opsim_dict, idx)
            self.modtran_visits.append(modtran_dict)

    def fillModtranDictionary(self, inputDict, idx):
        """Return a dictionary filled with all Modtran parameters
        for a given visit
        """
        if not self.atmos:
            raise ValueError('Atmosphere class not called')
        mdict = {}
        RA, DEC, mjd = (inputDict['fieldRA'],
                        inputDict['fieldDEC'],
                        inputDict['expMJD'])
        azimuth, z_angle = modtranTools.equatorial2local(RA, DEC,
                                                         mjd, unit='rad')
        mdict['ID'] = inputDict['obsHistID']
        mdict['ZANGLE'] = z_angle
        mdict['MODEL'] = self.atmos.model(idx)
        mdict['O3'] = self.atmos.ozone(mjd)
        mdict['H2O'] = self.atmos.vapor(mjd)
        mdict['VIS'] = self.getVIS(azimuth, z_angle, idx)
        mdict['ISEAS'] = self.atmos.iseas(idx)
        mdict['IVULC'] = self.atmos.ivulc(idx)
        mdict['IHAZE'] = self.atmos.ihaze()
        mdict['IVSA'] = self.atmos.ivsa()
        mdict['ZAER11'] = self.atmos.zaer11()
        mdict['ZAER12'] = self.atmos.zaer12()
        mdict['SCALE1'] = self.atmos.scale1(idx)
        mdict['ZAER21'] = self.atmos.zaer21()
        mdict['ZAER22'] = self.atmos.zaer22()
        mdict['SCALE2'] = self.atmos.scale2()
        return mdict

    def getVIS(self, azimuth, z_angle, idx):
        if not self.atmos:
            raise ValueError('Atmosphere class not called')
        if (z_angle*rad2dec) < 45.:
            sin2z = np.sin(2.*z_angle)
        else:
            sin2z = 1.
        vis0, visamp, visaz = (self.atmos.vis0[idx],
                               self.atmos.visamp[idx],
                               self.atmos.visaz[idx])
        vis = vis0 + visamp * np.sin(azimuth - visaz) * sin2z
        if vis < 4.:
            return 4.0
        else:
            return vis
