"""
modtranSequence.py

Class used to generate sequences of atmospheric parameters for input to MODTRAN.

Calling method :
--------------
>>> seq = AtmosphereSequence('Opsim3.61.02.sel.dat')
>>> seq.generateParameters(seed, 'atmos_database')

=>  generates parameters for this particular pointing with the provided seed
=>  saves the parameter list in $ATMOSPHERE_PARAMETERS_DIR/atmos_database.pck

>>> seq = AtmosphereSequence()
>>> seq.loadParameters('atmos_database.pck')
Parameters for # runs computed with seed = X
>>> seq.getAtmTrans('tmp', save=True)

=> computes the atmospheric extinction in $MODTRAN_DATADIR/tmp/
=> outputfiles in the form tmp.*
=> saves the final transmission data into $MODTRAN_DATADIR/tmp/tmp_final.plt

Description :
-----------
The AtmosphereSequence
- inputs a live opsim dictionary / list of dictionaries OR a filename,
- verifies the data contains at least the following keys :
'obsHistID', 'expMJD', 'fieldRA', 'fieldDEC',
- initializes a pointing sequence by defining the date boundaries,
- calls the Atmosphere class in the modtranAtmos.py and generates atmospheric
parameters for the given pointings,
- writes those parameters in a file / database (using the pickle module),
- loads parameters from database on demand,
- runs MODTRAN without the aerosols for each selected pointing (one at a time)
and returns output,
- computes the aerosol extinction spectrum for the same pointings,
- combines both spectra to output the final atmosphere transmission either
live or saved into a file.

To Lynne:
---------
NEED TO HAVE THE ENVIRONNEMENT VARIABLE DEFINED :
- LSST_POINTING_DIR targetting the path where the pointing files are.
- ATMOSPHERE_PARAMETERS_DIR targetting the parameter database

Need to check if MODTRAN output wavelengths are the same for version 4 and 5
=> cf. modwl.txt

"""

import numpy
import os
import pickle
from lsst.sims.atmosphere.transmission.modtranAtmos import Atmosphere
from lsst.sims.atmosphere.transmission.modtranCardsNoAer import ModtranCards
import lsst.sims.atmosphere.transmission.modtranTools as modtranTools

RAD2DEC = 180. / numpy.pi
DEG2RAD = numpy.pi / 180.

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
        """Instantiate an AtmosphereSequence."""
        # OpSim pointing file, if existing
        self.opsim_filename = opsim_filename
        # OpSim pointing dictionary, if existing
        self.opsim_data = opsim_data
        # List of OpSim visits for the sequence
        self.opsim_visits = []
        # List of MODTRAN dictionaries
        self.modtran_visits = []
        # List of aerosol parameters
        self.aerosol_visits = []
        # MODTRAN wavelengths array
        self.modtran_wl = numpy.array([])
        # Transmittance array
        self.transmittance = numpy.array([])

    def readOpsimData(self):
        """Read the pointing file and return the content as a list
        of dictionaries called opsim_visits. The dictionary keys are
        stored in class variable_opsim_keys.
        """
        if self.opsim_data:
            good, string = self.checkOpsimData()
            if good:
                self.opsim_visits = self.opsim_data
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
        else:
            raise Exception('No data specified')

    def checkOpsimData(self, otherdata=None):
        """Reads thru the data to see if there's any bug or missing key(s)."""
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
                statement = '{0} non dictionary element(s) and {1} uncomplete\
                    dictionary(ies) in data list'.format(nodict, baddict)
                return (False,  statement)

    def changeOpsimData(self, opsim_filename=None, opsim_data=None):
        """Input a new opsim data file or data list into the sequence."""
        self.opsim_filename = opsim_filename
        self.opsim_data = opsim_data

        self.opsim_visits = []

    def initPointingSequence(self):
        """Initialize input parameters for the atmospheric simulation."""
        if not self.opsim_visits:
            self.readOpsimData()
        # Number of visits
        self.npoints = len(self.opsim_visits)
        # Starting date
        self.mjds = self.opsim_visits[0]['expMJD']
        # Ending date
        self.mjde = self.opsim_visits[-1]['expMJD']

    def generateParameters(self, seed=_default_seed, output='atmos_db'):
        """Generate the atmospheric parameters over time.
        Returns a list of dictionaries containing the modtran
        information for each opsim visit (in the same order).
        """
        self.initPointingSequence()
        # Instantiate the Atmosphere class
        self.atmos = Atmosphere(
            self.mjds, self.mjde, self.npoints, seed)
        # Generate main atmosphere parameters sequence
        self.atmos.init_main_parameters()
        # Associate a value of these parameters for each pointing
        for opsim_dict in self.opsim_visits:
            # Get coordinates
            RA, DEC = (opsim_dict['fieldRA'], opsim_dict['fieldDEC'])
            # Get ID and date
            obsid, mjd = (opsim_dict['obsHistID'], opsim_dict['expMJD'])
            # Compute azimuth and elevation angle
            azimuth, z_angle = modtranTools.equatorial2local(RA, DEC,
                                                             mjd, unit='rad')
            # Get atmosphere parameters
            modtran_dict = self.fillModtranDictionary(mjd, obsid, z_angle)
            self.modtran_visits.append(modtran_dict)
            self.aerosol_visits.append(self.atmos.aerosols(mjd) + (z_angle,))
        if output:
            megatupl = (self.modtran_visits, self.aerosol_visits,)
            parmdir = os.getenv('ATMOSPHERE_PARAMETERS_DIR')
            outname = output + '.pck'
            parmpath = os.join.path(parmdir, outname)
            with open(parmpath, 'w') as parmf:
                pickle.dump(megatupl, parmf, seed)
        # Done

    def fillModtranDictionary(self, mjd, obsid, z_angle):
        """Return a dictionary filled with all Modtran parameters
        for a given visit."""
        if not self.atmos:
            raise ValueError('Atmosphere class not called')
        mdict = {}
        mdict['ID'] = obsid
        mdict['ZANGLE'] = z_angle
        mdict['MODEL'] = self.atmos.model(mjd)
        mdict['O3'] = self.atmos.ozone(mjd)
        mdict['H2O'] = self.atmos.wvapor(mjd)
        # put aerosols computation outsite of MODTRAN
        mdict['IHAZE'] = self.atmos.ihaze()
        return mdict

    def loadParameters(self, parmfile=''):
        """Load parameters for opsim visits to be computed with MODTRAN"""
        if not parmfile:
            raise IOError("You need to specify a parameter filename")
        parmdir = os.getenv('ATMOSPHERE_PARAMETERS_DIR')
        parmpath = os.join.path(parmdir, parmfile)
        # Read from file
        with open(parmpath, 'r') as parmf:
            data = pickle.load(parmf)
        # Dictionary list
        self.modtran_visits = data[0]
        # Tuple list
        self.aerosol_visits = data[1]
        # seed value
        nruns = len(self.modtran_visits)
        print 'Parameters for {1} runs computed with seed = {0}'.format(data[2],
                                                                        nruns)
        # Init transmission array
        self.initTransmissionArray(nruns)

    def getAtmTrans(self, outfile='tmp', save=True):
        """Go through the process of getting an atmospheric transmission"""
        # Set the name of MODTRAN output file
        self.outfilename = outfile
        # Loop over the MODTRAN runs
        for run in xrange(len(self.modtran_visits)):
            self.runModtran(run)
            # Read MODTRAN output
            modtrans = self.getModtranExtinction()
            # Compute aerosols transmission
            aertrans = self.getAerTransmittance(run)
            # Multiply both to get full transmission
            self.transmittance[run] = modtrans * aertrans

        if save:
            self.saveTrans()
        # Data saved into a file named '%s_final.plt' % self.outfilename

    def runModtran(self, run):
        """Call the modtranCards class in order to run MODTRAN for
        selected visits"""
        # Call for the MODTRAN card class
        modcard = ModtranCards()
        modcard.setDefaults()
        # Write the cards to the disk
        modcard.writeModtranCards(self.modtran_visits[run], self.outfilename)
        modcard.runModtran()

    def initTransmissionArray(self, nruns):
        """Initialize class attribute for atmospheric transmission"""
        if not self.modtran_wl:
            self.initModtranWavelengths()

        self.transmittance = numpy.zeros((nruns, len(self.modtran_wl)))

    def initModtranWavelengths(self):
        """Load in memory the sorted tabulated wavelengths in nanometers at which
        MODTRAN outputs data."""
        main_dir = os.getenv('ATMOSPHERE_TRANSMISSION_DIR')
        modtranwfile = os.path.join(main_dir, 'data/modtranwl.txt')
        self.modtran_wl = numpy.loadtxt(modtranwfile)

    def getModtranExtinction(self):
        """Read MODTRAN atmosphere extinction output file '.plt' and store the
        runs into a 2D array"""
        if not self.modtran_wl:
            self.initModtranWavelengths()

        modtranDataDir = os.getenv('MODTRAN_DATADIR')
        # MODTRAN transmission outputfile
        outputfile = '{0}/{1}.plt'.format(self.outfilename, self.outfilename)
        outputpath = os.path.join(modtranDataDir, outputfile)
        # Initialize array for transmittance
        trans = numpy.zeros(len(self.modtran_wl))
        with open(outputpath, 'r') as outf:
            # File starts with data - no header
            # I use negative indices because MODTRAN prints the wavelengths
            # in decreasing order
            idx = -1
            for line in outf:
                if line.starstwith('$'):
                    continue
                values = line.strip().split()
                trans[idx] = float(values[1])
                idx -= 1
                if abs(idx) > len(self.modtran_wl):
                    raise ValueError("Too many values to unpack from MODTRAN \
                        outputfile.")
        return trans
        # MODTRAN transmittance stored

    def saveTrans(self):
        """Save full extinction into an ascii file

        The structure is:
            wavelength transm[run1] transm[run2] transm[run3] etc.
        """
        modtranDataDir = os.getenv('MODTRAN_DATADIR')
        outputfile = '{0}/{1}_final.plt'.format(
            self.outfilename, self.outfilename)
        outputpath = os.path.join(modtranDataDir, outputfile)
        with open(outputpath, 'w') as transmf:
            transmf.write('$ FINAL ATMOSPHERE TRANSMISSION\n')
            for val in xrange(len(self.modtran_wl)):
                data = '\t'.join('{0:f}'.format(self.transmittance[run][val])
                                 for run in xrange(len(self.modtran_wl)))
                line = '{0}\t{1}\n'.format(self.modtran_wl[val], data)
                transmf.write(line)

    def getAerTransmittance(self, run):
        """Compute the atmospheric transmittance due to aerosols"""
        if not self.modtran_wl:
            self.initModtranWavelengths()

        # Compute Vandermonde matrix
        # degree is hardcoded but only 2nd degree polynomial is used
        vdm_wl = numpy.vander(self.modtran_wl, 3)
        # Get polynomial roots as well as zenith angle
        p0, p1, p2, z_ang = self.aerosol_visits[run]
        pfit = numpy.array([p0, p1, p2])
        # Reconstitute polynom
        polynom = numpy.dot(vdm_wl, pfit)
        # Retrieve airmass from zenith angle
        airmass = modtranTools.zenith2airmass(z_ang, site='lsst', unit='rad')
        return numpy.exp(-1.0 * airmass * polynom)
