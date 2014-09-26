####################################################################
#  Written by: Lynne Jones - UW - v1.1 5/21/10
#  Questions or comments, email : ljones@astro.washington.edu
#  $Id$
#  simSelfCalib.py
#  Python program to generate data for self-calibration trials
#
#    software requirements 
#         - numpy (version 1.2 or higher)
#         - mysqldb (optional - comment out *_db routines and import MySQLdb if not available)
#
#    read_parameters will read a parameter file, like simSelfCalib.input, 
#      and substitute defaults where inputs not specified
#    getVisits_db uses mysqldb to query opsim for the visits falling within 
#      ra/dec/time range specified (queries UW opsim db)
#    generateVisits creates a completely artifical sequence of visits to the
#      ra/dec range specified, with nEpoch return visits
#    addVisitErrors generates random zeropoint and color term errors for each visit
#    generateStars creates a completely artificial list of stars spread across 
#      ra/dec/mag/color range specified
#    getStars_db gets the stars from the database for the ra/dec/mag/color range specified
#    calcStarMags calculates the 'observed' magnitude for each star on each visit 
####################################################################



# assumptions: 
# information about observing schedule - visits and patches - will fit into memory
#  RA, Dec, time (or equiv to Nvisit), skyRot - note that RA and Dec are after taking                       
#   into account any dithering pattern.
# also zeropoint offset, zeropoint gradient, zeropoint gradient direction, and
# zeropoint changes due to color or a radial gradient in the color effect are stored
# in 'visits'.
#  These visits themselves are not currently broken down into patches, as the errors
# are smooth across the entire focal plane. Note that patches are solved independently, so
# this does not imply that the solver thinks the focal plane is smooth.
# The patches are used as estimates of how well the errors can be recovered though.
#  The stars are generated in blocks of the sky at one time.
#  NOTE - slight error with 'patch' calculation as patch can overlap different RA/DEc blocks,
# which is how the stars are generated ... so patchdata.dat can have 2 pieces of info for the same
# patch.

# general python modules
print 'importing general modules'
import numpy as np
import numpy.random as random
import time

# for illumination correction
import pyfits as pf
from scipy.interpolate import interp2d, interp1d

import healpy as hp

# importing modules for clouds
from scipy import interpolate
from scipy import stats
#print 'importing clouds'
import lsst.sims.selfcal.clouds.PowerSpectrum as pws
import lsst.sims.selfcal.clouds.Clouds as cld

import lsst.sims.maf.db as db



# If you have a problem with eups having two lsst.sims directories, please check that 
#  'lsst' is setup (should be, as part of your default setup, but please check)
#import PowerSpectrum as pws
#import Clouds as cld

#print 'importing useful I/O'
import lsst.sims.selfcal.analysis.useful_input as ui

import lsst.sims.selfcal.generation.FocalplaneThermalModel as fpm

# Need to fix this C-code to be able to import this
# import lsst.sims.selfcal.generation.ccs2amp_lib as cc

focalLength = 10300.0   # mm

rad2deg = 180.0/np.pi
deg2rad = np.pi/180.0
sqrt2 = np.sqrt(2)

# values for connecting to opsim database
# note this is likely a mysql db
#opsimhostname, opsimusername, opsimpassword, opsimdb, opsimdbtype = ui.setupConnEnv("LSST_OPSIM_")
#opsimtable = "output_opsim3_61"  

# values for connecting to calibration star database
# note this is a postgres db
#calsimhost, calsimusername, calsimpassword, calsimdb, calsimdbtype = ui.setupConnEnv("LSST_CALSIM_")
#calsimtable="msrgb_master"

# values for connecting to the local selfcal database
# note this is a postgres db
#selfcalhost, selfcalusername, selfcalpassword, selfcaldb, selfcaldbtype = ui.setupConnEnv("LSST_SELFCAL_")
#selfcaltable="selfcal"

# values for magnitude error estimation, fairly constant for LSST r band
#_systematicErr = 0.004
_systemM5 = 24.5

# default time start (for non-opsim visits)
_epochStart = 49353.
_epochStop = 49718.

verbose = True

# Utility function for distance on a sphere calculation.
def calcDist_cosines(RA1, Dec1, RA2, Dec2):
    """Calculates distance on a sphere using spherical law of cosines.
    Give this function RA/Dec values in radians. Returns angular distance(s), in radians.
    Note that since this is all numpy, you could input arrays of RA/Decs."""
    # This formula can have rounding errors for case where distances are small.
    # Oh, the joys of wikipedia - http://en.wikipedia.org/wiki/Great-circle_distance 
    # For the purposes of these calculations, this is probably accurate enough.
    D = np.sin(Dec2)*np.sin(Dec1) + np.cos(Dec1)*np.cos(Dec2)*np.cos(RA2-RA1)
    D = np.arccos(D)
    return D

def calcDist_haversine(RA1, Dec1, RA2, Dec2):
    """Calculates distance on a sphere using the haversine formula.
    Give this function RA/Dec values in radians. Returns angular distance(s), in radians.
    Note that since this is all numpy, you could input arrays of RA/Decs."""
    # This formula can have rounding errors for antipodal cases.
    D = (np.sin((Dec2-Dec1)/2))**2 + np.cos(Dec1)*np.cos(Dec2)*(np.sin((RA2-RA1)/2))**2
    D = np.sqrt(D)
    D = 2.0 * np.arcsin(D)
    return D

def calcDist_vincenty(RA1, Dec1, RA2, Dec2):
    """Calculates distance on a sphere using the Vincenty formula. 
    Give this function RA/Dec values in radians. Returns angular distance(s), in radians.
    Note that since this is all numpy, you could input arrays of RA/Decs."""
    # This one is supposed to be accurate for all distances, but is more complex. 
    D1 = (np.cos(Dec2)*np.sin(RA2-RA1))**2 + \
        (np.cos(Dec1)*np.sin(Dec2) - np.sin(Dec1)*np.cos(Dec2)*np.cos(RA2-RA1))**2
    D1 = np.sqrt(D1)
    D2 = (np.sin(Dec1)*np.sin(Dec2) + np.cos(Dec1)*np.cos(Dec2)*np.cos(RA2-RA1))
    D = np.arctan2(D1,D2)
    return D

def gnomonic_project_toxy(RA1, Dec1, RAcen, Deccen):
    """Calculate x/y projection of RA1/Dec1 in system with center at RAcen, Deccen.
    Input radians."""
    # also used in Global Telescope Network website
    cosc = np.sin(Deccen) * np.sin(Dec1) + np.cos(Deccen) * np.cos(Dec1) * np.cos(RA1-RAcen)
    x = np.cos(Dec1) * np.sin(RA1-RAcen) / cosc
    y = (np.cos(Deccen)*np.sin(Dec1) - np.sin(Deccen)*np.cos(Dec1)*np.cos(RA1-RAcen)) / cosc
    return x, y

def gnomonic_project_tosky(x, y, RAcen, Deccen):
    """Calculate RA/Dec on sky of object with x/y and RA/Cen of field of view.
    Returns Ra/Dec in radians."""
    denom = np.cos(Deccen) - y * np.sin(Deccen)
    RA = RAcen + np.arctan2(x, denom)
    Dec = np.arctan2(np.sin(Deccen) + y * np.cos(Deccen), np.sqrt(x*x + denom*denom))
    return RA, Dec
    
def simple_project(RA1, Dec1, RAcen, Deccen):
    """Calculate x/y projection of RA1/Dec from RAcen/Deccen, with simple projection.
    Don't use this - here for historical reasons. Input radians. """
    q = np.arctan(np.tan(Dec1)/np.cos(RA1-RAcen))
    x = np.cos(q) * np.tan(RA1-RAcen) / np.cos(q-RAcen)
    y = np.tan(q - Deccen)
    return x, y

def vignetFunc(x):
    """From VignettingFunc_v3.3.TXT.  r is in degrees, frac is fraction of rays which were not vignetted.
    returns the magnitudes of dimming caused by the vingetting relative to the center of the field"""

    if not hasattr(vignetFunc, 'r'):
        vignetFunc.r=np.array([0.000000,0.020000,0.040000,0.060000,0.080000,0.100000,0.120000,0.140000,
               0.160000,0.180000,0.200000,0.220000,0.240000,0.260000,0.280000,0.300000,
               0.320000,0.340000,0.360000,0.380000,0.400000,0.420000,0.440000,0.460000,
               0.480000,0.500000,0.520000,0.540000,0.560000,0.580000,0.600000,0.620000,
               0.640000,0.660000,0.680000,0.700000,0.720000,0.740000,0.760000,0.780000,
               0.800000,0.820000,0.840000,0.860000,0.880000,0.900000,0.920000,0.940000,
               0.960000,0.980000,1.000000,1.020000,1.040000,1.060000,1.080000,1.100000,
               1.120000,1.140000,1.160000,1.180000,1.200000,1.220000,1.240000,1.260000,
               1.280000,1.300000,1.320000,1.340000,1.360000,1.380000,1.400000,1.420000,
               1.440000,1.460000,1.480000,1.500000,1.520000,1.540000,1.560000,1.580000,
               1.600000,1.620000,1.640000,1.660000,1.680000,1.700000,1.720000,1.740000,
               1.760000,1.780000,1.800000,1.820000,1.840000,1.860000,1.880000,1.900000,
               1.920000,1.940000,1.960000,1.980000,2.000000])
        vignetFunc.frac=np.array([0.623885,0.623885,0.623885,0.623885,0.623885,0.623885,0.623885,0.623885,
                  0.623885,0.623885,0.623885,0.623885,0.623885,0.623885,0.623885,0.623885,
                  0.623885,0.623885,0.623885,0.623885,0.623885,0.623885,0.623885,0.623885,
                  0.623885,0.623885,0.623885,0.623885,0.623885,0.623885,0.623885,0.623885,
                  0.623885,0.623822,0.623759,0.623632,0.623442,0.623316,0.623000,0.622494,
                  0.622367,0.621861,0.621671,0.621292,0.621039,0.620659,0.620216,0.619963,
                  0.619394,0.619204,0.618635,0.618319,0.618066,0.617496,0.617117,0.616737,
                  0.616168,0.615535,0.614840,0.614207,0.613385,0.612436,0.611614,0.610602,
                  0.609716,0.608957,0.608071,0.606996,0.605668,0.604972,0.603770,0.602758,
                  0.601177,0.599595,0.598140,0.595673,0.594282,0.592447,0.590613,0.588526,
                  0.586312,0.584287,0.582137,0.580113,0.578025,0.576064,0.573344,0.570561,
                  0.563730,0.545322,0.521412,0.492757,0.460624,0.429565,0.404959,0.383073,
                  0.356190,0.317161,0.279777,0.241824,0.201974])
        vignetFunc.frac=vignetFunc.frac/vignetFunc.frac[0]

    result = -2.5*np.log10(np.interp(x,vignetFunc.r,vignetFunc.frac))
    return result

def keplerStats(nstars):
    """ Use distributions from McQuillan et al 2011 to assign variability magnitude to each star.
    Note, this does not assign things properly based on color, etc. """
    #r_var in mmag
    r_var = np.array([  1.72023666e-01,   4.34372473e-01,   4.68398291e-01,
         4.93441626e-01,   5.14732312e-01,   5.33943909e-01,
         5.52017102e-01,   5.67281763e-01,   5.83079492e-01,
         5.98434571e-01,   6.12578979e-01,   6.26858191e-01,
         6.40081335e-01,   6.53214771e-01,   6.66383496e-01,
         6.79295044e-01,   6.92332715e-01,   7.05689294e-01,
         7.19002673e-01,   7.32294849e-01,   7.45732666e-01,
         7.59157483e-01,   7.71666126e-01,   7.83999268e-01,
         7.97270813e-01,   8.09827271e-01,   8.23078759e-01,
         8.35341333e-01,   8.48190100e-01,   8.61307434e-01,
         8.74334473e-01,   8.87676636e-01,   9.00621155e-01,
         9.13179846e-01,   9.25284668e-01,   9.37659424e-01,
         9.50595105e-01,   9.62950037e-01,   9.75506982e-01,
         9.87822754e-01,   9.99954407e-01,   1.01253831e+00,
         1.02530010e+00,   1.03845610e+00,   1.05139209e+00,
         1.06403748e+00,   1.07686282e+00,   1.08966384e+00,
         1.10221550e+00,   1.11597666e+00,   1.12970129e+00,
         1.14343447e+00,   1.15815095e+00,   1.17344751e+00,
         1.18902659e+00,   1.20497607e+00,   1.22192629e+00,
         1.23907109e+00,   1.25790491e+00,   1.27619675e+00,
         1.29680090e+00,   1.31895554e+00,   1.34181226e+00,
         1.36509580e+00,   1.39152437e+00,   1.42032690e+00,
         1.45045315e+00,   1.48488755e+00,   1.52108020e+00,
         1.55989465e+00,   1.60238257e+00,   1.64887463e+00,
         1.69844570e+00,   1.75147295e+00,   1.81425400e+00,
         1.87909448e+00,   1.94843804e+00,   2.03531321e+00,
         2.12714468e+00,   2.22709292e+00,   2.32715478e+00,
         2.43474858e+00,   2.54564619e+00,   2.67783809e+00,
         2.84044277e+00,   3.06426636e+00,   3.36279756e+00,
         3.73507793e+00,   4.22938799e+00,   4.77742559e+00,
         5.42256055e+00,   6.15611787e+00,   7.05875137e+00,
         8.11473389e+00,   9.28510273e+00,   1.08929326e+01,
         1.28745992e+01,   1.55171701e+01,   1.94736664e+01,
         2.74170559e+01,   4.73471844e+02])
    percent_cumulative = np.array([ 0.  ,  0.01,  0.02,  0.03,  0.04,  0.05,  0.06,  0.07,  0.08,
        0.09,  0.1 ,  0.11,  0.12,  0.13,  0.14,  0.15,  0.16,  0.17,
        0.18,  0.19,  0.2 ,  0.21,  0.22,  0.23,  0.24,  0.25,  0.26,
        0.27,  0.28,  0.29,  0.3 ,  0.31,  0.32,  0.33,  0.34,  0.35,
        0.36,  0.37,  0.38,  0.39,  0.4 ,  0.41,  0.42,  0.43,  0.44,
        0.45,  0.46,  0.47,  0.48,  0.49,  0.5 ,  0.51,  0.52,  0.53,
        0.54,  0.55,  0.56,  0.57,  0.58,  0.59,  0.6 ,  0.61,  0.62,
        0.63,  0.64,  0.65,  0.66,  0.67,  0.68,  0.69,  0.7 ,  0.71,
        0.72,  0.73,  0.74,  0.75,  0.76,  0.77,  0.78,  0.79,  0.8 ,
        0.81,  0.82,  0.83,  0.84,  0.85,  0.86,  0.87,  0.88,  0.89,
        0.9 ,  0.91,  0.92,  0.93,  0.94,  0.95,  0.96,  0.97,  0.98,
        0.99,  1.  ])
    percent_draws = np.random.rand(nstars)
    r_var_out = np.interp(percent_draws,  percent_cumulative,  r_var)
    return r_var_out

def variType(r_var):
    """Assign a variability type to each star.  The types are 1 = Gaussian, 2 = rotational.
    For now, I'll leave flaring stars out of it, but that's another possibility"""
    #from the paper, 16% of all dwarfs have rotational variability, and I think most of those are on the high variability tail
    split_pt = 1.6 #mmag, an arbitrary pull between low and high variability
    high_var = np.ravel(np.where(r_var > split_pt)) #define high variability as 2x solar
    low_var = np.ravel(np.where(r_var < split_pt)) #complement
    type=np.zeros(np.size(r_var))+1
    #make an array where 90% are 2s and the rest are 1s
    temp=np.random.rand(np.size(high_var))
    type[high_var[np.where(temp > 0.5)]]=2 #set 50% of the high variability stars to be periodic
    temp=np.random.rand(np.size(low_var))
    type[low_var[np.where(temp > 0.98)]]=2 #set 2% of the low variability stars to be periodic
    return type

def variOffset(r_var, type):
    """ given the magnitude of variabilty r_var (5-95 percentile range) and type of variability
    1= Gaussian, 2 = rotational/periodic, generate random offsets.  This will be too variable
    for periodic variability and observations that are closely spaced in time."""
    sigma = 3.3*r_var
    dmag = np.zeros(np.size(r_var)+0.)
    g = np.ravel(np.where(type == 1))
    dmag[g]=np.random.randn(np.size(g))*sigma[g]
    per = np.ravel(np.where(type ==2))
    dmag[per]=0.5*r_var[per]*np.sin(2.*np.pi*np.random.rand(np.size(per)))
    return dmag/1000. #convert to mags

def illumcorrFunc(x, y):
    val = np.zeros((len(x)))

    if not hasattr(illumcorrFunc, 'illumcorrInterp'):
        return val
    else:
        for i in range(len(val)):
            r = np.sqrt(x[i]**2 + y[i]**2)
            val[i] = illumcorrFunc.illumcorrInterp(r)

            if val[i] <= 0:
                val[i] = 0
            else:
                val[i] = -2.5*np.log10(val[i])
        
        return val

def process_illumcorr(filename):
    """ The illumination correction file is a text file containing the illumination correction
    as a ratio (not a magnitude).  The illumination correction is assumed to be only a function of r [0,1]
    """

    try:
        icdat = np.loadtxt(filename)
    except IOError:
        raise IOError

    print "illumcorrInterp"

    interpFunc = interp1d(icdat[:,0], icdat[:,1], bounds_error=False, fill_value=0)

    illumcorrFunc.illumcorrInterp = interpFunc
    
    print "done with process_illumcorr"

    return
    

def wrapRA(ra):
    """Wraps RA into 0-360 degrees."""
    ra = ra % 360.0
    return ra

def capDec(dec):
    """Terminates declination at +/- 90 degrees."""
    dec = np.where(dec>90, 90, dec)
    dec = np.where(dec<-90, -90, dec)
    return dec
        
def readParameters(parameter_filename):
    """ Read input parameters for simSelfCalib.py output generation.
    
    Parameter file inputs should be keyword = value. See simSelfCalib.input for example. """
    # These are default values, but can be reset through parameter file.
    # These describe the telescope setup - field of view and # of patches per field of view. 
    radius_fov = 1.8
    nPatch = 16
    # These describe the number of visits and the number of stars. 
    nEpoch = 10
    nStarTot = 1000000
    id_start = 0
    # These describe the footprint of the survey on the sky, and the properties of the stars.
    raMin = 0
    raMax = 360
    decMin = -80.
    decMax = 0.0
    colmax = 3.5
    colmin = -0.5
    magmin = 17
    magmax = 21
    random_seed = None
    errorDist = 'Gaussian'
    systematicErr = 0.004
    fluxNoise = 0.01
    fluxFrac = 0
    # These describe the errors to add to the simulation model. 
    mag_rand_err = 0.01
    calcerror = True
    zp_var_max = 1.0
    zp_grad_max = 0.3
    colterm_max = 0.2
    colterm_rad_max = 0.1
    color_obs_noise = 0.05
    color_correction_error=0.10
    filter_jitter = 0.005
    gainvar=0.00
    exptvar = 0.00
    nside = 8
    healshift = 2.3
    use_cloudsimage = False
    cloud_scale = 0.034
    sinvarx_mag = 0.0
    sinvary_mag = 0.0
    sinvarx_scale = 1.0
    sinvary_scale = 1.0
    sinx_phase = 0.0
    siny_phase = 0.0
    sinx_angle = 0.0
    flat_sinvarx_mag = 0.00
    flat_sinvary_mag = 0.00
    flat_sinvarx_scale = 1.
    flat_sinvary_scale = 1.
    flat_sinx_angle = 0.0
    phase_persist = 1
    cloud_mag = 0.10
    cloud_sinvar_scale = 5.
    rr_fraction = 0
    rr_amplitude = 1.
    illum_corr_filename = None
    fpModel_filename = None
    # These describe how the field location changes between visits (visit=epoch). 
    dith_raOff_frac = 0.5
    dith_decOff_frac = 0.5
    dith_skyRot_min = -90
    dith_skyRot_max = 90
    # These describe an alternate method of setting the field locations. 
    use_opsim = False
    use_opsimdither = False
    opsimfilter = "r"
    tstart = 0
    tstop = 360
    # These describe an alternate method of setting the star locations.
    use_calsim = False
    calsimtable = "bhb"
    # use a local database for output
    dbOut = False
    #
    kepler_variability = False
    #
    # Parameters associated with focalplane temperature variation - opsimfilter must be 'y'
    # for these to have any effect
    #
    fpModel1_filename = None
    fpModel2_filename = None
    fpTempTimescale = 0.2 # days
    ymagTempColorSlope = 0.001 * (29.3 - 25.8)/(6.0 * 5.0) # d(y4mag)/d(temp)d(g-i color) mag per degK per mag
    ymagTempSlope = 0.001 * 26.2 / 5.0 # mag per degK at g-i=0
    #
    # These describe the default data output file names.
    starobs_filename = "star_obs.dat"
    print_subpatch = False
    master_filename = "master_cal.dat"
    stardata_filename = "stardata.dat"
    patch_filename = "patchdata.dat"
    visit_filename = "visit.dat"
    illumcorr_filename = None
 
    
    # Open the parameter file and read inputs (to override values above, if given). 
    if (parameter_filename == "NULL"):
        print "No parameter file given, going with the defaults"
        return (raMin, raMax, decMin, decMax, use_opsim, use_opsimdither, opsimfilter,
                use_calsim, calsimtable, dbOut,
                radius_fov, nPatch, nEpoch, tstart, tstop,
                mag_rand_err, calcerror, zp_var_max, zp_grad_max, 
                colterm_max, colterm_rad_max, color_obs_noise,color_correction_error,
                filter_jitter,gainvar,exptvar,nside,healshift,use_cloudsimage,cloud_scale,sinvarx_mag, sinvary_mag, 
                sinvarx_scale, sinvary_scale, sinx_phase, siny_phase, sinx_angle,
                flat_sinvarx_mag, flat_sinvary_mag, flat_sinvarx_scale, 
                flat_sinvary_scale, flat_sinx_angle, phase_persist,
                cloud_mag, cloud_sinvar_scale, rr_fraction, rr_amplitude, 
                dith_raOff_frac, dith_decOff_frac, dith_skyRot_min, dith_skyRot_max,
                magmin, magmax, nStarTot, colmin, colmax, random_seed, errorDist, systematicErr, fluxNoise, fluxFrac,  
                starobs_filename, print_subpatch, master_filename, stardata_filename, patch_filename, visit_filename, 
                illumcorr_filename, fpModel1_filename, fpModel2_filename, fpTempTimescale, ymagTempColorSlope,
                ymagTempSlope)
    try:
        parfile = open(parameter_filename, 'r')
    except IOError:
        print "Could not open parameter file %s, going with defaults" %(parameter_filename)
        return (raMin, raMax, decMin, decMax, use_opsim, use_opsimdither, opsimfilter, 
                use_calsim, calsimtable, dbOut,
                radius_fov, nPatch, nEpoch, tstart, tstop,
                mag_rand_err, calcerror,zp_var_max, zp_grad_max, colterm_max, colterm_rad_max, color_obs_noise,color_correction_error,
                filter_jitter,gainvar,exptvar,nside,healshift,use_cloudsimage, cloud_scale,sinvarx_mag,
                sinvary_mag, sinvarx_scale, sinvary_scale, sinx_phase, siny_phase, sinx_angle,
                flat_sinvarx_mag, flat_sinvary_mag, flat_sinvarx_scale, flat_sinvary_scale, flat_sinx_angle, phase_persist,
                cloud_mag, cloud_sinvar_scale, rr_fraction, rr_amplitude,
                dith_raOff_frac, dith_decOff_frac, dith_skyRot_min, dith_skyRot_max,
                magmin, magmax, nStarTot, colmin, colmax, random_seed, errorDist, systematicErr, fluxNoise, fluxFrac, 
                starobs_filename, print_subpatch, master_filename, stardata_filename, patch_filename, visit_filename,
                illumcorr_filename, fpModel1_filename, fpModel2_filename, fpTempTimescale, ymagTempColorSlope,
                ymagTempSlope)
    # Read the parameter file inputs.
    for lines in parfile: 
        if lines.startswith("#"):
            continue
        keyword, junk, value = lines.split()
        if keyword == "nStarTot":
            nStarTot = int(value)
        elif keyword == "id_start":
            id_start = int(value)
        elif keyword == "magmin":
            magmin = float(value)
        elif keyword == "magmax":
            magmax = float(value)
        elif keyword == "colmax":
            colmax = float(value)
        elif keyword == "colmin":
            colmin = float(value)
        elif keyword == "random_seed":
            random_seed = int(value)
        elif keyword == "errorDist":
            errorDist = value
        elif keyword == "systematicErr":
            systematicErr = float(value)
        elif keyword == "fluxNoise":
            fluxNoise = float(value)
        elif keyword == "fluxFrac":
            fluxFrac = float(value)
        elif keyword == "mag_rand_err":
            mag_rand_err = float(value)
        elif keyword == "calcerror":
            calcerror = value
      # Try to account for a few different spellings. 
            if (calcerror=="TRUE") | (calcerror=="true") | (calcerror=="True") | (calcerror=="yes"):
                calcerror = True
            elif (calcerror=="FALSE") | (calcerror=="false") | (calcerror=="False") | (calcerror=="no"):
                calcerror = False
            else :
                print "Sorry I don't recognize that value for calcerror - please try True or False"
                exit()      
        elif keyword == "zp_var_max":
            zp_var_max = float(value)
        elif keyword == "zp_grad_max":
            zp_grad_max = float(value)
        elif keyword == "colterm_max":
            colterm_max = float(value)
        elif keyword == "colterm_rad_max":
            colterm_rad_max = float(value)
        elif keyword == "color_obs_noise":
            color_obs_noise = float(value)
        elif keyword == "color_correction_error":
            color_correction_error = float(value)
        elif keyword == "filter_jitter":
            filter_jitter = float(value)
        elif keyword == "gainvar":
            gainvar = float(value)
        elif keyword == "exptvar":
            exptvar = float(value)
        elif keyword == "nside":
            nside = int(value)
        elif keyword == "healshift":
            healshift = float(value)
        elif keyword == "use_cloudsimage":
            use_cloudsimage  = value
      # Try to account for a few different spellings. 
            if (use_cloudsimage=="TRUE") | (use_cloudsimage=="true") | (use_cloudsimage=="True") | (use_cloudsimage=="yes"):
                use_cloudsimage = True
            elif (use_cloudsimage=="FALSE") | (use_cloudsimage=="false") | (use_cloudsimage=="False") | (use_cloudsimage=="no"):
                use_cloudsimage = False
            else :
                print "Sorry I don't recognize that value for use_cloudsimage - please try True or False"
                exit()                  
        elif keyword == "cloud_scale":
            cloud_scale = float(value)
        elif keyword == "sinvarx_mag":
            sinvarx_mag = float(value)
        elif keyword == "sinvary_mag":
            sinvary_mag = float(value)
        elif keyword == "sinvarx_scale":
            sinvarx_scale = float(value)
        elif keyword == "sinvary_scale":
            sinvary_scale = float(value)
        elif keyword == "sinx_phase":
            sinx_phase = float(value)
        elif keyword == "siny_phase":
            siny_phase = float(value)
        elif keyword == "sinx_angle":
            sinx_angle = float(value)
        elif keyword == "flat_sinvarx_mag":
            flat_sinvarx_mag = float(value)
        elif keyword == "flat_sinvary_mag":
            flat_sinvary_mag = float(value)
        elif keyword == "flat_sinvarx_scale":
            flat_sinvarx_scale = float(value)
        elif keyword == "flat_sinvary_scale":
            flat_sinvary_scale = float(value)
        elif keyword == "flat_sinx_phase":
            flat_sinx_phase = float(value)
        elif keyword == "flat_sinx_angle":
            flat_sinx_angle = float(value)
        elif keyword == "phase_persist":
            phase_persist = float(value)
        elif keyword == "cloud_mag":
            cloud_mag= float(value)
        elif keyword == "cloud_sinvar_scale":
            cloud_sinvar_scale = float(value)
        elif keyword == "rr_fraction":
            rr_fraction = float(value)
        elif keyword == "rr_amplitude":
            rr_amplitude = float(value)
        elif keyword == "kepler_variablity":
            kepler_variability = value
            if (kepler_variability=="TRUE") | (kepler_variability=="true") | (kepler_variability=="True") | (kepler_variability=="yes"):
                kepler_variability = True
            elif (kepler_variability=="FALSE") | (kepler_variability=="false") | (kepler_variability=="False") | (kepler_variability=="no"):
                kepler_variability = False
            else :
                print "Sorry I don't recognize that value for kepler_variability - please try True or False"
                exit()            
        elif keyword == "raMin":
            raMin = float(value)
        elif keyword == "raMax":
            raMax = float(value)
        elif keyword == "decMin":
            decMin = float(value)
        elif keyword == "decMax":
            decMax = float(value)
        elif keyword == "use_opsim" : 
            use_opsim = value
            # Try to account for a few different spellings. 
            if (use_opsim=="TRUE") | (use_opsim=="true") | (use_opsim=="True") | (use_opsim=="yes"):
                use_opsim = True
            elif (use_opsim=="FALSE") | (use_opsim=="false") | (use_opsim=="False") | (use_opsim=="no"):
                use_opsim = False
            else :
                print "Sorry I don't recognize that value for use_opsim - please try True or False"
                exit()
        elif keyword == "use_opsimdither":
            use_opsimdither = value
            if ((use_opsimdither=="TRUE") | (use_opsimdither=="true") | (use_opsimdither=="True")
                | (use_opsimdither=="yes")):
                use_opsimdither = True
            elif ((use_opsimdither=="FALSE") | (use_opsimdither=="false") |
                  (use_opsimdither=="False") | (use_opsimdither=="no")):
                use_opsimdither = False
        elif keyword == "opsimfilter":
            opsimfilter = value
            if ((opsimfilter=="u") | (opsimfilter=="g") | (opsimfilter=="r") |
                (opsimfilter=="i") | (opsimfilter=="z") | (opsimfilter=="y")):
                pass
            else:
                print "Sorry, I don't recognize that filter choice for opsim: please try u g r i z or y"
                exit()
        elif keyword == "use_calsim" : 
            use_calsim = value
            # Try to account for a few different spellings. 
            if (use_calsim=="TRUE") | (use_calsim=="true") | (use_calsim=="True") | (use_calsim=="yes"):
                use_calsim = True
            elif ((use_calsim=="FALSE") | (use_calsim=="false") | (use_calsim=="False")
                  | (use_calsim=="no")):
                use_calsim = False
            else :
                print "Sorry I don't recognize that value for use_calsim - please try True or False"
                exit()
        elif keyword == "calsimtable":
            calsimtable = value
        elif keyword == "dbOut":
            dbOut = value
      # Try to account for a few different spellings. 
            if (dbOut=="TRUE") | (dbOut=="true") | (dbOut=="True") | (dbOut=="yes"):
                dbOut = True
            elif (dbOut=="FALSE") | (dbOut=="false") | (dbOut=="False") | (dbOut=="no"):
                dbOut = False
            else :
                print "Sorry I don't recognize that value for dbOut - please try True or False"
                exit()      
        elif keyword == "tstart":
            tstart = float(value)
        elif keyword == "tstop":
            tstop = float(value)
        elif keyword == "nEpoch":
            nEpoch = int(value)
        elif keyword == "radius_fov":
            radius_fov = float(value)
        elif keyword == "nPatch":
            nPatch = int(value)
            pside = int(np.sqrt(nPatch))
            print " nPatch implies %d x %d grid on field of view" %(pside, pside)
            if (np.abs(pside*pside  - nPatch) > 0.1):
                print " Please make input nPatch a square number (i.e. 16, as 16=4x4) to remove ambiguity"
                exit()
        elif keyword == "dith_raOff_frac":
            dith_raOff_frac = float(value)
        elif keyword == "dith_decOff_frac":
            dith_decOff_frac = float(value)
        elif keyword == "dith_skyRot_min":
            dith_skyRot_min = float(value)
        elif keyword == "dith_skyRot_max":
            dith_skyRot_max = float(value)
        elif keyword == "starobs_filename":
            starobs_filename = value
        elif keyword == "print_subpatch":
            print_subpatch = value
            if (print_subpatch=="TRUE") | (print_subpatch=="true") | (print_subpatch=="True") | (print_subpatch=="yes"):
                print_subpatch = True
            elif (print_subpatch=="FALSE") | (print_subpatch=="false") | (print_subpatch=="False") | (print_subpatch=="no"):
                print_subpatch = False
            else :
                print "Sorry I don't recognize that value for print_subpatch - please try True or False"
                exit()      
        elif keyword == "master_filename":
            master_filename = value
        elif keyword == "stardata_filename":
            stardata_filename = value
        elif keyword == "patch_filename":
            patch_filename = value
        elif keyword == "visit_filename":
            visit_filename = value
        elif keyword == "illumcorr_filename":
            illumcorr_filename = value
        elif keyword == "fpModel1_filename":
            fpModel1_filename = value
        elif keyword == "fpModel2_filename":
            fpModel2_filename = value
        elif keyword == "ymagTempSlope":
            ymagTempSlope = float(value)
        elif keyword == "ymagTempColorSlope":
            ymagTempColorSlope = float(value)
        elif keyword == "fpTempTimescale":
            fpTempTimescale = float(value)
        else:
            print value, " Sorry I don't recognize that keyword. Please check the code or example input file."
            exit()
    # Done reading information from parameter file.
    parfile.close()
    # Just do a couple of checks.
    # Check that RA/Dec values are within bounds. 
    # Handle dec limits in a rather limiting way, since this is hard to decipher. 
    decMax = capDec(decMax)
    decMin = capDec(decMin)
    # Handle ra limits by wrapping them into 0-360 degrees. 
    # Note that if raMin > raMax, we will then have to assume that this means the pointing was 
    # intended to wrap through RA=0.  
    raMin = wrapRA(raMin)
    raMax = wrapRA(raMax)
    # Check that time to 'stop' survey is larger than starting time.
    if tstart > tstop:
        junk = tstop
        tstop = tstart
        tstart = junk
    if use_opsim==True:
        print "Using opsim to generate visits, with filter choice of %s" %(opsimfilter)
        if use_opsimdither==True:
            print "  and using opsim's hex dither pattern." 
        if (tstart == tstop):
            print "You want to use the opsim to generate visits, but tstart = tstop. Please fix."
            exit()
    if use_opsim==False:
        print "Using the regular grid to generate %d visits to each pointing in the grid" %(nEpoch)
        if nEpoch == 0:
            print "You want to generate visits without opsim, but the number of repeat visits is 0. "
            print "Please fix this problem in the input file." 
            exit()
    if use_calsim==True:
        print "Using calibration database %s.%s on %s to generate stars." %(calsimdb,
                                                                            calsimtable, calsimhost)
    else:
        print "Using a flat random distribution of %d stars. " %(nStarTot)

    # If an illumination correction file has been specified, process it
    if illumcorr_filename:
        try:
            process_illumcorr(illumcorr_filename)
        except IOError:
            print "Could not process illumcorr_file %s" % illumcorr_filename
            exit()

    # Return values read from input file. 
    return (raMin, raMax, decMin, decMax, use_opsim, use_opsimdither, opsimfilter,
            use_calsim, calsimtable, dbOut,
            radius_fov, nPatch, nEpoch, tstart, tstop,
            mag_rand_err, calcerror,zp_var_max, zp_grad_max, colterm_max, colterm_rad_max, color_obs_noise,color_correction_error,
            filter_jitter,gainvar,exptvar,nside,healshift,use_cloudsimage,cloud_scale,sinvarx_mag,
            sinvary_mag, sinvarx_scale, sinvary_scale, sinx_phase, siny_phase, sinx_angle,
            flat_sinvarx_mag, flat_sinvary_mag, flat_sinvarx_scale, flat_sinvary_scale, flat_sinx_angle, phase_persist,
            cloud_mag, cloud_sinvar_scale, rr_fraction, rr_amplitude,kepler_variability, dith_raOff_frac, dith_decOff_frac, dith_skyRot_min, dith_skyRot_max,
            magmin, magmax, nStarTot, id_start, colmin, colmax, random_seed, errorDist, systematicErr, fluxNoise, fluxFrac, 
            starobs_filename, print_subpatch, master_filename, stardata_filename, patch_filename, visit_filename,
            illumcorr_filename, fpModel1_filename, fpModel2_filename, fpTempTimescale, ymagTempColorSlope,
            ymagTempSlope)


def getVisits_db(ra_min, ra_max, dec_min, dec_max, time_start, time_stop, opsimfilter, 
                 rad_fov_deg, use_opsimdither, dith_raOff_frac, dith_decOff_frac, random_seed):
    """Get visits by querying opsim database. 
    
    Return all visits falling within ra/dec boundary, within given time.
    Sky rotation comes from opsim, but dithering of field of view comes from 
      dith_raOff_frac, dith_decOff_frac.
      Time_start/stop is in days (since start of survey), ra/dec in degrees. """
    # Connect to the database.
    conn, cursor = ui.sqlConnect(hostname=opsimhostname, username= opsimusername, passwdname = opsimpassword, 
                                 dbname=opsimdb, dbtype=opsimdbtype)
    # Opsim uses radians for field locations, so convert and add buffer (fudge factor) for radius of
    #  field of view because want to include fields which may fall partly within this boundary. 
    # This buffer here includes a sqrt2 because I only want to include fields which have the inscribed
    # square overlapping the desired field of view (rather than the tangent point of the circle). 
    decmin = dec_min - rad_fov_deg/sqrt2 
    decmax = dec_max + rad_fov_deg/sqrt2
    decmin = capDec(decmin)
    decmax = capDec(decmax)
    decmin = decmin*deg2rad
    decmax = decmax*deg2rad
    maxdra_of_rad_fov = np.abs(rad_fov_deg / np.cos(max(np.abs(decmin), np.abs(decmax))) / sqrt2)
    ramin = (ra_min - maxdra_of_rad_fov) 
    ramax = (ra_max + maxdra_of_rad_fov)    
    ramin = wrapRA(ramin)
    ramax = wrapRA(ramax)
    # Check to see if addition of rad_fov_deg pushed ramin/max past each other, like if
    # input ramin/max was 0/360, then with rad_fov_deg, would now be 351/8 or so.
    print "RA limits, before and after wrapping: ", ramin, ramax, ra_min, ra_max
    if (ramin>=ramax) & (ra_min<=ra_max):
        ramin=0
        ramax=360
    ramin = ramin*deg2rad
    ramax = ramax*deg2rad
    print "Querying for fields from %s with ra between %f and %f, dec between %f and %f" \
              %(opsimtable, ramin*rad2deg, ramax*rad2deg, decmin*rad2deg, decmax*rad2deg)
    if use_opsimdither:
        sqlquery = 'select DISTINCT ON (expmjd) hexdithra, hexdithdec, rotSkyPos, night, expmjd, "5sigma_modified" from %s' %(opsimtable)
    else:
        sqlquery = 'select DISTINCT ON (expmjd) fieldra, fielddec, rotSkyPos, night, expmjd, "5sigma_modified" from %s' %(opsimtable)
    # declination is in order from -90 to 90
    sqlquery = sqlquery + " where (fielddec between %f and %f)" %(decmin, decmax)
    # ra might wrap around RA = 0 though
    if (ramin<ramax):
        sqlquery = sqlquery + " and (fieldRA between %f and %f)" % (ramin, ramax)
    else: 
        sqlquery = sqlquery + " and ((fieldRA>%f) or (fieldRA<%f))" %(ramin, ramax)
    sqlquery = sqlquery + " and (night between %f and %f)" % (time_start, time_stop)
#    sqlquery = sqlquery + " and filter='%s' group by expmjd order by min(expmjd)" %(opsimfilter)
    sqlquery = sqlquery + " and filter='%s'  " %(opsimfilter)
    print "# ", sqlquery
    # Get Query results.
    sqlresults = ui.sqlQuery(cursor, sqlquery)
    print "Fetched %d potential fields from database" %(len(sqlresults))
    ui.sqlEndConnect(conn, cursor)
    # Assign query results to appropriate fields in dictionary.
    visits = {}
    keys = ('ra', 'dec', 'skyrot', 'night', 'time', '5sigma_modified')
    arraylen = len(sqlresults)
    if (arraylen==0):
        raise Exception("Found no observations matching these constraints.")
    dtypes = ['float' ,'float' ,'float' ,'float', 'float','float' ,'int' ]
    for i in np.arange(np.size(keys)):
        visits[keys[i]] = np.empty((arraylen,), dtype=dtypes[i])
    i = 0
    for result in sqlresults: 
        j = 0  
        for key in keys:
            visits[key][i] = result[j]
            j = j+1
        i = i+1
    # Convert RA/dec back to degrees.
    visits['ra'] = visits['ra'] * rad2deg
    visits['dec'] = visits['dec'] * rad2deg    
    # Add dithering to RA/Dec of each visit (unless use_opsimdither is true). 
    # The dithering is a flat random distribution, within the limits specified by dith_*Off_frac.
    # Generate array of RA and Dec offsets (random.rand returns 'nvisits' random numbers between 0-1).
    if not (use_opsimdither):
        random.seed(random_seed)
        RAoffsets = rad_fov_deg * dith_raOff_frac * random.rand(arraylen)
        DecOffsets = rad_fov_deg * dith_decOff_frac * random.rand(arraylen)
        # Apply these RA and Dec offsets to the RA/Dec central pointings. 
        visits['dec'] = visits['dec'] + DecOffsets
        visits['ra'] = visits['ra'] + RAoffsets / np.cos(visits['dec']*deg2rad)
    # Check for any wrap-around/under 0/360 RA or dec after dither
    visits['dec'] = capDec(visits['dec'])
    visits['ra'] = wrapRA(visits['ra'])
    # Sky rotation also in radians, so convert also.
    visits['skyrot'] = visits['skyrot'] * rad2deg
    # Should also cull visits which actually don't fall within ra/dec limits + radius_fov (because
    # could have added too many by using too-large a fudge factor earlier) or were dithered out of limits
    ramin = wrapRA(ra_min)
    ramax = wrapRA(ra_max)
    decmin = capDec(dec_min)
    decmax = capDec(dec_max)
    #for i in range(arraylen):
    if (ra_min < ra_max):
        condition = ((visits['ra'] >= ramin - rad_fov_deg/np.cos(visits['dec']*deg2rad)/sqrt2) &
                     (visits['ra'] <= ramax + rad_fov_deg/np.cos(visits['dec']*deg2rad)/sqrt2) &
                     (visits['dec'] >= decmin - rad_fov_deg/sqrt2) &
                     (visits['dec'] <= decmax + rad_fov_deg/sqrt2))
    else: # ra_min > ra_max, wrapping through ra=0
        condition = (((visits['ra'] >= ramin - rad_fov_deg/np.cos(visits['dec']*deg2rad)/sqrt2) |
                      (visits['ra'] <= ramax + rad_fov_deg/np.cos(visits['dec']*deg2rad)/sqrt2)) &
                     (visits['dec'] >= decmin - rad_fov_deg/sqrt2) &
                     (visits['dec'] <= decmax + rad_fov_deg/sqrt2))
    visits['ra'] = visits['ra'][condition]
    visits['dec'] = visits['dec'][condition]
    visits['time'] = visits['time'][condition]
    visits['skyrot'] = visits['skyrot'][condition]
    visits['night'] = visits['night'][condition]
    visits['5sigma_modified']=visits['5sigma_modified'][condition]
    nvisits = len(visits['ra'])
    visits['id'] = np.arange(0, nvisits, 1, dtype='int')
    print "Got %d visits from the opsim over the desired RA/Dec and time range" %(nvisits)
    return visits

def generateVisits(ra_min, ra_max,  dec_min, dec_max, nepochs, rad_fov_deg, 
                   dith_raOff_frac, dith_decOff_frac, dith_skyRot_min, dith_skyRot_max, random_seed):
    """Generate visits by simply tiling observation area to number of epochs/revisits, using given FOV.
    
    RA, Dec, FOV are input in degrees. The RA/decOff_frac are fractions of the total field of view 
    which the fields are dithered over. The skyRot_min/max gives the minimum and maximum camera+telescope
    rotation angle (skyRot = sky rotation angle .. =angle between North and the +Y axis in the camera),
    in degrees.  The tiled area returned may be *larger* than the area requested, as it is fully covered,
    but no field centers (undithered) will be returned that are beyond the requested boundaries. 
    The separation of each center of the field of view is not 2*radius, but (radius/sqrt(2) * 2).
    Returns arrays of RA/Dec/skyRotation for each visit. """
    # dith_skyrot_max = Maximum sky rotation angle (angle between north and +Y axis), in degrees,
    # while dith_skyrot_min = minimum sky rotator angle. Values of (-90, 90) are consistent with LSST.
    # The final size of this array would be *something* like nepochs*d(RA)*d(dec)/rad_fov_deg^2.`
    visits = {}
    visits['ra'] = np.empty((0), dtype='float')
    visits['dec'] = np.empty((0), dtype='float')
    # Set up regular grid through RA/Dec.
    ra_grid = np.empty((0), dtype='float')
    dec_grid = np.empty((0), dtype='float')
    # Calculate how many dec bins we'll actually need.
    decmin = capDec(dec_min)
    decmax = capDec(dec_max)
    nbins_dec = int(np.ceil((decmax - decmin) / (2.0*rad_fov_deg/sqrt2)))
    ramin = wrapRA(ra_min)
    ramax = wrapRA(ra_max)
    if ramax == ramin:
        ramax = ramin + 360.0
    # Step through dec grid, build RA steps. 
    for bin_dec in range(0, nbins_dec, 1):
        dec = decmin + (bin_dec+0.5)*(rad_fov_deg*2.0/sqrt2)
        if (dec > decmax) & (bin_dec>0):
            continue
        # Set the location of the RA center of field of view, including the factor of cos(dec)
        #  so that the fields are separated by the radius of the field of view. 
        rad_fov_deg_here = rad_fov_deg / np.cos(dec*deg2rad)
        if (ramin < ramax):
            ra_row = np.arange(ramin, ramax+rad_fov_deg_here/sqrt2, 2.0*rad_fov_deg_here/sqrt2)
            # in case ramin+rad_fov_deg_here is greater than ramax already ... just add one fov at this dec
        else: # have to wrap through RA=0
            ra_maxTemp = ramax + 360.0
            # allocate using a grid with values larger than 360
            ra_row = np.arange(ramin, ra_maxTemp+rad_fov_deg_here/sqrt2,
                              2.0*rad_fov_deg_here/sqrt2)
        # then check didn't go past 360 degrees
        ra_row = wrapRA(ra_row)
        if (len(ra_row)==0):
            ra_row = np.zeros((1), dtype='float')
            ra_row[0] = (ramin + ramax) / 2.0
        # Calculate the declination for this row of RA.
        dec_row = np.zeros(len(ra_row), dtype='float')
        dec_row = dec_row + dec
        # And add these onto the grid (the reason for the resize to 0 at the start).
        ra_grid = np.append(ra_grid, ra_row)
        dec_grid = np.append(dec_grid, dec_row)
    # Replicate grid for each visit.
    visits['ra'] = ra_grid
    visits['dec'] = dec_grid
    for visit in range (1, nepochs, 1):  
        visits['ra'] = np.append(visits['ra'], ra_grid)
        visits['dec'] = np.append(visits['dec'], dec_grid)    
    # And add dithering to RA/Dec of each visit. 
    # The dithering is a flat random distribution, within the limits specified by dith_*Off_frac.
    nvisits = len(visits['ra'])
    # Generate array of RA and Dec offsets (random.rand returns 'nvisits' random numbers between 0-1).
    random.seed(random_seed)
    RAoffsets = rad_fov_deg * dith_raOff_frac * random.rand(nvisits)
    DecOffsets = rad_fov_deg * dith_decOff_frac * random.rand(nvisits)
    # Apply these RA and Dec offsets to the RA/Dec central pointings. 
    visits['dec'] = visits['dec'] + DecOffsets
    visits['ra'] = visits['ra'] + RAoffsets / np.cos(visits['dec']*deg2rad)
    # Check the RA and Dec limits
    visits['ra'] = wrapRA(visits['ra'])
    visits['dec'] = capDec(visits['dec'])
    # Time is a little meaningless in this case, but set up some values.
    timestep = (_epochStop - _epochStart) / len(visits['ra'])
    visits['time'] = np.arange(_epochStart, _epochStop, timestep)
    # Generate a flat distribution of sky rotation angles.
    visits['skyrot'] = (dith_skyRot_max - dith_skyRot_min) * random.rand(nvisits) + dith_skyRot_min
    # Should also cull visits which actually don't fall within ra/dec limits + radius_fov
    # were dithered out of limits.
    arraylen = len(visits['ra'])
    print "Looking for visits between RA %f and %f, Dec %f and %f" %(ramin, ramax, decmin, decmax)
    #for i in range(arraylen):        
    if (ramin < ramax):
        condition = ((visits['ra'] >= ramin - rad_fov_deg/np.cos(visits['dec']*deg2rad)/sqrt2) &
                     (visits['ra'] <= ramax + rad_fov_deg/np.cos(visits['dec']*deg2rad)/sqrt2) &
                     (visits['dec'] >= decmin - rad_fov_deg/sqrt2) &
                     (visits['dec'] <= decmax + rad_fov_deg/sqrt2))
    else: # ra_min > ra_max, wrapping through ra=0
        condition = (((visits['ra'] >= ramin - rad_fov_deg/np.cos(visits['dec']*deg2rad)/sqrt2) |
                      (visits['ra'] <= ramax + rad_fov_deg/np.cos(visits['dec']*deg2rad)/sqrt2)) &
                     (visits['dec'] >= decmin - rad_fov_deg/sqrt2) &
                     (visits['dec'] <= decmax + rad_fov_deg/sqrt2) )
    visits['ra'] = visits['ra'][condition]
    visits['dec'] = visits['dec'][condition]
    visits['time'] = visits['time'][condition]
    visits['skyrot'] = visits['skyrot'][condition]
    visits['id'] = visits['id'][condition]
    nvisits = len(visits['ra'])
    visits['id'] = np.arange(0, nvisits, 1, dtype='int')
    #need to add in 'night' and '5sigma_modified' parameters to keep things compatible
    visits['night']=np.round(visits['time'])
    visits['night']=visits['night']-np.min(visits['night'])
    visits['5sigma_modified']=visits['time']*0+24.6826 #just use an average r-band depth
    print "Generated %d visits from the opsim over the desired RA/Dec and time range" %(nvisits)
    return visits

def createVisitErrors(visits, zp_max, zp_grad_max, colterm_max, colterm_rad_max, filter_jitter, visit_filename,
                      random_seed, phase_persist, gainvar):
    """Add zeropoint and color term errors to each patch within each visit. 
    
    zp=overall offset, zp_grad=linear slope in zp, colterm=linear coefficient (mag/mag), 
       colterm_radgrad=radial dependence for colterm (mag/mag/deg). 
    Returns ra/dec for each visit with dither, zeropoint/zeropointgrad, color/colorgrad. """
    random.seed(random_seed)
    nvisits = len(visits['ra'])
    # Generate flat distribution of random zeropoint offsets between zero and zp_max. 
    #  note that random.rand(nvisits) returns a n array of values between 0-1 of length nvisits.
    visits['zpOff'] = zp_max*random.rand(nvisits)
    # Generate flat distribution of random zeropoint gradients. (0-1)
    visits['zpGrad'] = zp_grad_max * random.rand(nvisits)
    # Generate random angle for zeropoint gradient.
    visits ['zpGradDir']= 360.0*random.rand(nvisits)
    # Generate overall color term (flat distribution). (0-1)
    visits['colorterm'] = colterm_max * random.rand(nvisits)
    # Generate radial color term slope (flat random distribution). (0-1)
    visits['colortermRad'] = colterm_rad_max/2.*np.ones(nvisits) #taking a factor of two here to convert it to a slope.
    # Generate a fixed patern for flat-field for each night.
    visits['flatPhase'] = np.zeros(nvisits)
    # make the days blocks of "phase_persist" length, so the pattern can repeat for multiple days.  
    days=np.ceil(np.around(visits['night'])/phase_persist)
    udays = np.unique(days)
    uphases = random.rand(len(udays))*2.*np.pi
    namp=7200 #30x240
    if gainvar != 0:
        visits['ampvar']=np.zeros((nvisits,namp))
    for i in np.arange(len(udays)):
        good=np.where(days == udays[i])
        visits['flatPhase'][good]=uphases[i]
        if gainvar != 0:
            visits['ampvar'][good]=visits['ampvar'][good]+np.random.randn(namp)
    # assume there has been a filter change if the time between observations is above 135 s 
    tsort = np.argsort(visits['time'])
    tgap = np.roll(visits['time'][tsort],-1)-visits['time'][tsort]
    tgap[0] = 0
    visits['filterPosition'] = np.zeros(len(tgap))
    change_time = 120. + 15. #seconds, from the Science Book, leaves out slew time, but that should be small
    change_positions = np.where(tgap > change_time/(60.*60.*24.))
    visits['filterPosition'][change_positions]= 1. * random.randn(len(np.ravel(change_positions)))*filter_jitter
    visits['angles']=np.zeros(len(tgap))
    visits['angles'][change_positions]=random.rand(len(np.ravel(change_positions)))*2*np.pi
    #loop through
    visits['filterXoff'] = visits['filterPosition']*0
    visits['filterYoff'] = visits['filterPosition']*0
    if filter_jitter != 0:
        for i in range(np.min(np.where( visits['filterPosition'] > 0)), len(visits['filterPosition'])):
            if visits['filterPosition'][i] == 0 :
                visits['filterPosition'][i] = visits['filterPosition'][i-1]
                visits['angles'][i]=visits['angles'][i-1]
    #re-arange so that things are back in the right order
        visits['filterPosition'][tsort] = visits['filterPosition']
        visits['angles'][tsort]=visits['angles']
    #angles=random.rand(len(visits['filterPosition']))*2*np.pi
        visits['filterXoff'] = visits['filterPosition'] * np.cos(visits['angles'])
        visits['filterYoff'] = visits['filterPosition'] * np.sin(visits['angles'])
    print_visitfile(vfile=None, visit_output=visit_filename, visits=visits)
    return visits

def generateStars(sfile, id_start, ra_min, ra_max, dec_min, dec_max, mag_min, mag_max, nstars, 
                  col_min=-0.5, col_max=3.5, random_seed=None):
    """Generate a random (flat) mag, color, ra, and dec distribution of stars. 

    Returns stars within ra_min/max, dec_min/max space ONLY - does not add stars to cover full fov.
    ra and dec are input in degrees."""
    # Check range of dec
    decmin = capDec(dec_min)
    decmax = capDec(dec_max)    
    # Check range of RA.
    ramin = wrapRA(ra_min)
    ramax = wrapRA(ra_max)
    if (ramin >= ramax):
        ramax = ramax + 360.0
    # Set random RA within ramin/ramax range, and sort so stars are in RA order.
    stars = {}
    random.seed(random_seed)
    stars['ra']= (ramax - ramin) * random.rand(nstars) + ramin
    stars['ra'].sort()
    # Set random dec within decmin/decmax range, with flat distribution on sphere.
    # use sin instead of cos because matches with natural range of n sin
    vmax = (np.sin(decmax*deg2rad)+1.0)/2.0
    vmin = (np.sin(decmin*deg2rad)+1.0)/2.0
    v = (vmax-vmin)*random.rand(nstars) + vmin
    stars['dec'] = np.arcsin(2*v-1)
    stars['dec'] = stars['dec']*rad2deg
    # Wrap ra, in case they were split through RA=0.
    stars['ra'] = wrapRA(stars['ra'])
    # set random magnitude distribution (flat) and color distribution
    stars['rmagtrue'] = (mag_max - mag_min) * random.rand(nstars) + mag_min
    stars['color'] = (col_max-col_min) * random.rand(nstars) + col_min
    # give each star a unique identification number
    stars['id'] = np.arange(id_start, id_start+nstars, 1, dtype=int)
    stars['dbid'] = stars['id']
    # write star information out to file, for future reference (also in master file)
    sfile = print_starstruth(sfile, stars=stars)
    # return values
    id_next = stars['id'][len(stars['id'])-1] + 1
    return stars, nstars, id_next

def getStars_db(sfile, cursor, calsimtable, id_start, ra_min, ra_max, dec_min, dec_max,
                mag_min, mag_max, nstars, col_min=-0.5, col_max=3.5, filter='r'):
    """Get star information from star database (pass cursor to db connection). 

    Writes stars true magnitudes and colors to a flat file, when pulled from database.
    Returns stars within ra_min/max, dec_min/max ONLY, ra and dec are in degrees. """
    #XXX time check
    #t1=time.time()
    requested_mag=filter+'mag'
    #print 'TIME CHECK starting db query %f'%t1
    # Verify inputs. 
    ramin = wrapRA(ra_min)
    ramax = wrapRA(ra_max)
    decmin = capDec(dec_min)
    decmax = capDec(dec_max)
    magmax = mag_max
    # Construct query.
    queryobj = 'count(rmag)'
    queryobj = 'simobjid, ra, decl, %s, umag, gmag, rmag, imag, zmag, ymag'%requested_mag #bring this up to just skip the count check
    querycount = 1
    loopcount = 0
    while querycount<2:
        sqlquery = "select %s from %s" %(queryobj, calsimtable)
        sqlquery = sqlquery + " where (%smag between %f and %f)" %(filter,mag_min, magmax)
        sqlquery = sqlquery + " and (decl between %f and %f)" %(decmin, decmax)
        if ra_min < ra_max:
            sqlquery = sqlquery + " and ((ra >= %f) and (ra < %f))" %(ramin, ramax)
        else:
            sqlquery = sqlquery + " and ((ra >= %f) or (ra < %f))" %(ramin, ramax)
        if querycount > 0:
            sqlquery = sqlquery + " order by rmag limit %i" %(nstars)
            #sqlquery = sqlquery + " order by ra"
        #if verbose:
        #    print "# ", sqlquery
        cursor.execute(sqlquery)
        sqlresults = cursor.fetchall()
        loopcount += 1
        # If we've already queried for simobjid/ra/.., quit.
        if querycount == 1:
            querycount = 2
        # Otherwise, see if count is within limits. 
        if querycount == 0:
            counts = int(sqlresults[0][0])
            querycount = 1
            ##Removing this chunk because it can blow up if the query is in a spot where the MW is cut out.
#            if (counts>nstars/1.5) | (loopcount>5):
#                querycount = 1
#                queryobj = 'simobjid, ra, decl, rmag, gmag, imag'
#            else:
#                factor = (magmax - mag_min)*(nstars/counts)
#                magmax = mag_min + factor*(magmax - mag_min)
    # Got query results. 
    #t2=time.time()
    #print 'TIME CHECK got db results in %f min'%((t2-t1)/60.)
    stars = {}
#    keys = ('dbid', 'ra', 'dec', 'rmagtrue', 'gmag', 'imag')
    keys = ('dbid', 'ra', 'dec', 'rmagtrue', 'umag', 'gmag', 'rmag', 'imag', 'zmag', 'ymag')
    arraylen = len(sqlresults)
    # Bail out if no stars were returned from database.
    if arraylen==0:
        nstars = 0
        for key in keys:
            stars[key] = None  ### 
        print 'no stars'
        return stars, nstars, id_start
    # Otherwise assign query results to dictionary of arrays.
    for key in keys:
        stars[key] = np.empty((arraylen,), dtype='float')
    i = 0
    for result in sqlresults:
        j = 0
        for key in keys:
            stars[key][i] = result[j]
            j = j+1
        i = i+1
    # Trim stars beyond the requested color range.
    print 'I wanted nstars = ', nstars  ####XXX
    print 'I got nstars =', len(stars['ra']) ###XXX
    stars['color'] = stars['gmag'] - stars['imag']
    condition = (stars['color'] > col_min) & (stars['color'] < col_max)
    keys = keys + ('color',)    
    for key in keys:
        stars[key] = stars[key][condition]
    nstars = len(stars['ra'])
    print 'after color-cut, nstars =', nstars #XXX
    # give each star a unique identification number
    stars['id'] = np.arange(id_start, id_start+nstars, 1, dtype=int)
    id_next = stars['id'][len(stars['id'])-1] + 1
    # Write stars information out to file
    sfile = print_starstruth(sfile, stars=stars)
    # Return stars.
    #t2=time.time()
    #print 'TIME CHECK retuning db results in %f min'%((t2-t1)/60.)
    return stars, nstars, id_next

def generateCloudImage(corrFunc,rad_fov_deg, sampling=256):
    windowsize=np.sqrt(2.*(2.*rad_fov_deg)**2.)
    #ps=pws.PowerSpectrum(windowsize,sampling)
    #ps.ComputeStructureFunction()
    #ps.ComputeCorrelationFunction()
    c=cld.Clouds(windowsize,sampling)
    #c.DirectClouds(ps.getCorrel2D())
    c.DirectClouds(corrFunc)
    x,y=np.mgrid[0:sampling,0:sampling]
    return x,y,c.clouds


def generateStarMags(visits, ra_min, ra_max, dec_min, dec_max, mag_min, mag_max,
                     col_min, col_max, color_obs_noise, color_correction_error,cloud_mag, cloud_sinvar_scale, \
                     rr_fraction, rr_amplitude,variability, rad_fov_deg, fieldscale, \
                     ofile, print_subpatch, mfile, sfile, pfile, fpTempModel1, fpTempModel2,
                     mag_rand_error, gainvar, exptvar, nside, healshift, calcerror,
                     npatch=16, use_calsim=False, use_cloudsimage=False,  cloud_scale=0.034,
                     cursor=None, calsimtable=None, nstars=None,
                     id_start=0, random_seed=None, fluxNoise=0.01, fluxFrac=0, filter='r', errorDist='Gaussian', systematicErr=0.004):
    """Given visits, generates magnitudes for each observation of a star within ra/dec min/max.
    
    Pulls information for stars from a database using getStars_db or from random star distribution using
    generateStars. Calculates their magnitudes in each visit, stacking visits to the same area of sky
    in order to reduce db queries and only generate stars once.
    Call this routine once for a particular ra/dec region. 
    Writes out values for star's magnitudes as they are calculated.  """
    # Find visits that overlap this block of sky.
    ramin = wrapRA(ra_min)
    ramax = wrapRA(ra_max)
    decmin = capDec(dec_min)
    decmax = capDec(dec_max)
    if (ramin < ramax):
        condition = ( (visits['ra'] >= ramin - rad_fov_deg/np.cos(visits['dec']*deg2rad)) &
                      (visits['ra'] <= ramax + rad_fov_deg/np.cos(visits['dec']*deg2rad)))
    else:
        condition = ( (visits['ra'] >= ramin - rad_fov_deg/np.cos(visits['dec']*deg2rad)) |
                      (visits['ra'] <= ramax + rad_fov_deg/np.cos(visits['dec']*deg2rad)))
    condition = condition & (visits['dec'] > decmin - rad_fov_deg) & (visits['dec'] < decmax + rad_fov_deg)
    blockvisits = {}
    for key in visits.keys():
        blockvisits[key] = visits[key][condition]
    if len(blockvisits['ra']) == 0:
        return id_start
    # Get stars that fall into this block of sky. 
    if use_calsim:
        stars, nstars_found, id_next = getStars_db(sfile, cursor, calsimtable, id_start,
                                                   ra_min, ra_max, dec_min, dec_max,
                                                   mag_min, mag_max, nstars, col_min, col_max, filter=filter)
    else:        
        stars, nstars_made, id_next = generateStars(sfile, id_start, ra_min, ra_max, dec_min, dec_max,
                                                    mag_min, mag_max, nstars, col_min, col_max,
                                                    random_seed)
    #if we hit a patch where there are no stars, just bail out
    if (blockvisits['ra'] == None) or (stars['ra'] == None):
        return id_next
    if verbose:
        print "In RA %f:%f, Dec %f:%f, found %d compatible visits, with %d possible star matches." \
              %(ramin, ramax,
                decmin, decmax,
                len(blockvisits['ra']), len(stars['ra']))
    #Generate flux standard stars
    if fluxFrac != 0:
        starsCAL = {}
        if fluxFrac < 1:
            nCalStars = len(stars['ra'])*fluxFrac 
            if nCalStars < 1:
                if np.random.rand(1) < nCalStars:
                    nCalStars = 1
                else:
                    nCalStars = 0            
        else:
            nCalStars = np.round(fluxFrac)
        print 'Making %f flux stars in block'%nCalStars
        if nCalStars != 0:
            condition = np.round(np.random.rand(nCalStars)*(nCalStars-1))
            condition = condition.astype(int)
            starsCAL['id'] = stars['id'][condition]
            starsCAL['fullpatch'] = starsCAL['id']*0 #just set all the flux standards to patch zero
            starsCAL['rmagobs'] = stars['rmagtrue'][condition] +  fluxNoise*np.random.randn(nCalStars)
            starsCAL['magerr'] = stars['rmagtrue'][condition]*0+fluxNoise
            ofile = print_obsfile(ofile, stars=starsCAL)
            pfile = print_patchfile(pfile, patch={'visitid':0, 'fullpatch':0, 'subpatch':0,
                                                  'count':nCalStars, 'ra':0, 'dec':0, \
                                                  'rmagtrue':0, 'color':0, \
                                                  'dmag':fluxNoise, 'dmag_rms':fluxNoise, \
                                                  'dmag_rand':0, 'dmag_rand_rms':0, \
                                                  'dmag_zp':0, 'dmag_zp_rms':0, \
                                                  'dmag_color':0, 'dmag_color_rms':0})
            starsCAL['dbid']=starsCAL['subpatch']=starsCAL['dmag_var']= starsCAL['dmag_snr']= \
                              starsCAL['dmag_zp']= starsCAL['dmag_zp_dist']= \
                              starsCAL['dmag_color']= starsCAL['dmag_color_dist']= starsCAL['dmag_sin']= \
                              starsCAL['dmag_flat_sin']= starsCAL['dmag_illumcorr']=starsCAL['dmag_cloud']=starsCAL['dmag_cloud_image']= \
                              starsCAL['dmag_rr']=starsCAL['rmagtrue']= starsCAL['color']= starsCAL['ra']= \
                              starsCAL['dec']= \
                              starsCAL['X']= starsCAL['Y']= starsCAL['dm5']=starsCAL['dmag_kep']= \
                              starsCAL['dmag_gainvar']=starsCAL['dmag_exptvar']=  starsCAL['rmagobs']*0
            starsCAL['night'] = 0
            mfile = print_masterfile(mfile, visitid = blockvisits['id'][0], raCen=0, decCen=0,
                                     stars=starsCAL)                                    
    #generate the power spectrum for the clouds.  
    if use_cloudsimage:
        sampling=256#512#256#128
        print "using cloud sampling = %i"%sampling
        windowsize=np.sqrt(2.*(2.*rad_fov_deg)**2.)
        ps=pws.PowerSpectrum(windowsize,sampling)
        ps.ComputeStructureFunction()
        ps.ComputeCorrelationFunction()
        c=cld.Clouds(windowsize,sampling)
        powerspec=ps.getCorrel2D()
    # clip the variability dictionary down to only include stars in the visit
    if variability:
        sub_variability = variability.copy()
        stars_in_visit=np.in1d(variability['id'], stars['id'])
        for key in sub_variability.keys():
            sub_variability[key]=sub_variability[key][stars_in_visit]
    # Go through each visit in block & calculate what stars fall within fov and where (rough x/y).
    for visitN in range(0, len(blockvisits['ra'])):
        raCen =  blockvisits['ra'][visitN]
        decCen = blockvisits['dec'][visitN]
        m5 = blockvisits['5sigma_modified'][visitN]
        # Calculate the angular distances of each star from the center of the field of view.
        starDist = rad2deg * calcDist_cosines(raCen*deg2rad, decCen*deg2rad, 
                                              stars['ra']*deg2rad, stars['dec']*deg2rad)
        condition = (starDist < rad_fov_deg)
        # Pull out the stars within this field of view - these are still arrays.
        starsFOV = {}
        for key in stars.keys():
            starsFOV[key] = stars[key][condition]
        # Asign the stars to a HEALpixel
        #HPdec = (decCen+90.)*deg2rad
        #HPra = raCen*deg2rad
        HPdec = (starsFOV['dec']+90.)*deg2rad
        HPra = starsFOV['ra']*deg2rad
        #HP1 = hp.ang2pix(nside,HPdec,HPra)
        #starsFOV['HPid1'] = starsFOV['id']*0+HP1
        starsFOV['HPid1'] = hp.ang2pix(nside,HPdec,HPra)
        #hps, hp_weights = hp.get_neighbours(nside,HPdec,phi=HPra)
        hps = hp.get_neighbours(nside,HPdec,phi=HPra)[0]
        hp_keys = ['HPid2','HPid3','HPid4', 'HPid5']
        for i in np.arange(np.size(hp_keys)):
            starsFOV[hp_keys[i]] = hps[i]
        #print 'XXX, size of HP arrays', np.size(starsFOV['ra']), np.size(starsFOV['HPid1']), np.size(starsFOV['HPid2'])
#        for key in hp_keys:
#            starsFOV[key] = starsFOV['id']*0-1
#        weight_cutoff = 0.001
#        good = np.ravel(np.where( (hp_weights > weight_cutoff) & (hps != HP1)))
#        for i in np.arange(np.size(good)):
#            starsFOV[hp_keys[i]] = starsFOV['id']*0+hps[good[i]]

        # Calculate the star's offset from center of this field of view (in square).
        starsFOV['X'], starsFOV['Y'] = gnomonic_project_toxy(starsFOV['ra']*deg2rad,
                                                             starsFOV['dec']*deg2rad,
                                                             raCen*deg2rad, decCen*deg2rad)        
        # And include rotation angle of visit.
        # skyRot angle is the angle between Y towards North .. so use -skyRotangle
        sin_skyrot = np.sin(blockvisits['skyrot'][visitN]*deg2rad)
        cos_skyrot = np.cos(blockvisits['skyrot'][visitN]*deg2rad)
        x = np.copy(starsFOV['X']) # make a copy of the x coord
        starsFOV['X'] = cos_skyrot*x + sin_skyrot*starsFOV['Y']
        starsFOV['Y'] = -1.*sin_skyrot*x + cos_skyrot*starsFOV['Y']
        # Get the raft, ccd, and ccd pixel coords for the stars
        starsFOV['ccd_Name'], starsFOV['ccd_X'], starsFOV['ccd_Y'] = fp_map_xy(starsFOV['X'], starsFOV['Y'])

        if opsimfilter == 'y':
            # Get the temperature at that focalplane position
            a1 = 0.5 * (1.0 + np.sin(blockvisits['time'][visitN]/fpTempTimescale))
            a2 = 1. - a1
            fpTemp1 = fp_map_temp(fpTempModel1, starsFOV['ccd_Name'],starsFOV['ccd_X'], starsFOV['ccd_Y'], blockvisits['time'][visitN])
            fpTemp2 = fp_map_temp(fpTempModel2, starsFOV['ccd_Name'],starsFOV['ccd_X'], starsFOV['ccd_Y'], blockvisits['time'][visitN])
            starsFOV['dtemp'] = a1*fpTemp1 + a2*fpTemp2
            starsFOV['dmag_temp'] = starsFOV['dtemp']*(ymagTempSlope + ymagTempColorSlope*starsFOV['color'])
        else:
            starsFOV['dtemp'] = np.zeros((len(starsFOV['X'])))
            starsFOV['dmag_temp'] = np.zeros((len(starsFOV['X'])))

        # starRadius is in degrees
        starRadius = np.sqrt(starsFOV['X']**2+starsFOV['Y']**2)/fieldscale*rad_fov_deg
        #apply vignetting
        starsFOV['dm5']=vignetFunc(starRadius)
        m5=m5-starsFOV['dm5']
        #Record which night we are on
        starsFOV['night']=blockvisits['night'][visitN]
        # Then calculate color term errors - must first calculate radial color offset distance. 
        colterm_dist = np.sqrt((starsFOV['X']+blockvisits['filterXoff'][visitN])**2 + \
                              (starsFOV['Y']+blockvisits['filterYoff'][visitN])**2) / fieldscale
        #make a "measured color" term      
        starsFOV['dmag_color_dist'] = colterm_dist*blockvisits['colortermRad'][visitN] * \
                                      (starsFOV['color'])
        #if we want to remove the color gradient noise now:
        if color_obs_noise > 0:
            ctd_unshifted = np.sqrt(starsFOV['X']**2 + starsFOV['Y']**2) / fieldscale
            starsFOV['measured_color']=starsFOV['color']+color_obs_noise*random.randn(len(starsFOV['color']))          
            starsFOV['dmag_color_dist'] =starsFOV['dmag_color_dist'] -  ctd_unshifted*blockvisits['colortermRad'][visitN] * \
                                      (starsFOV['measured_color'])+\
                                      random.randn(len(starsFOV['color']))*color_correction_error*starsFOV['dmag_color_dist'] 
        #add the color term noise
        starsFOV['dmag_color'] = blockvisits['colorterm'][visitN] * (starsFOV['color'])
        # Next calculate gray zero point offset - gradient plus whole fov.
        # Need to calculate distance along zp gradient direction for zp grad error. 
        zpGrad_dist = (starsFOV['X'] * np.cos(blockvisits['zpGradDir'][visitN]*deg2rad) +
                       starsFOV['Y'] * np.sin(blockvisits['zpGradDir'][visitN]*deg2rad)) / fieldscale
        starsFOV['dmag_zp_dist'] = zpGrad_dist * blockvisits['zpGrad'][visitN]
        starsFOV['dmag_zp'] = blockvisits['zpOff'][visitN] * np.ones(len(starsFOV['X']), dtype='float')
        # Calculate the sinusoidal variation x' = x*cos(phi) + y*sin(phi)
        xp=starsFOV['X']*np.cos(sinx_angle)+starsFOV['Y']*np.sin(sinx_angle)
        starsFOV['dmag_sin'] = sinvarx_mag*np.sin(np.pi*xp/fieldscale*sinvarx_scale + sinx_phase) +\
                               sinvary_mag*np.sin(np.pi*starsFOV['Y']/fieldscale*sinvary_scale + siny_phase)     
        #XXX -- so blockvisits['flatPhase'][visitN] should work
        xp=starsFOV['X']*np.cos(flat_sinx_angle)+starsFOV['Y']*np.sin(flat_sinx_angle)
        starsFOV['dmag_flat_sin'] = flat_sinvarx_mag*np.sin(np.pi*xp/fieldscale*flat_sinvarx_scale + \
                                                           blockvisits['flatPhase'][visitN] ) +\
                               flat_sinvary_mag*np.sin(np.pi*starsFOV['Y']/fieldscale*flat_sinvary_scale + \
                                                      blockvisits['flatPhase'][visitN])
        # interpolate illumination correction
        starsFOV['dmag_illumcorr'] = illumcorrFunc(starsFOV['X']/fieldscale, starsFOV['Y']/fieldscale)
        
        # Put in some structure for the clouds
        total_cloud_mag = cloud_mag*random.rand(1)*blockvisits['zpOff'][visitN]
        cloud_phase = random.rand(1)*2.*np.pi
        starsFOV['dmag_cloud'] = total_cloud_mag/2.*np.sin(np.pi*starsFOV['X']/fieldscale*cloud_sinvar_scale+cloud_phase)+\
                                 total_cloud_mag/2.*np.sin(np.pi*starsFOV['Y']/fieldscale*cloud_sinvar_scale+cloud_phase)
        #Make more realistic clouds
        if use_cloudsimage:
            t1=time.time()
            cloudx,cloudy,cloudImage=generateCloudImage(powerspec,rad_fov_deg, sampling=sampling)
            #cloudx,cloudy,cloudImage=generateCloudImage(powerspec, rad_fov_deg, sampling=256)
            t2=time.time()
            cloudx=cloudx/(np.max(cloudx)/2.)-1
            cloudy=cloudy/(np.max(cloudy)/2.)-1
            cloudxy=np.zeros( (np.size(cloudx),2))
            cloudxy[:,0]=np.ravel(cloudx)
            cloudxy[:,1]=np.ravel(cloudy)  #there clearly has to be a better way to do that...
            cloudImage = cloudImage/np.std(cloudImage)*cloud_scale*blockvisits['zpOff'][visitN]
            #cloudImage=cloudImage-np.min(cloudImage)
            #scale the cloudImage to be between 0.5 and 1.5, then multiply by the mean zeropoint offset.
            #No idea if this is a reasonable range.  
            #cloudImage=(cloudImage/np.max(cloudImage)+0.5)*blockvisits['zpOff'][visitN]-blockvisits['zpOff'][visitN]#XXX??!?!?!?? no idea how to scale this image properly subtract off the zpOff again so that the base zeropint doesn't get appl
            tempxy=np.zeros( (np.size(starsFOV['X']),2))
            tempxy[:,0]=starsFOV['X']/fieldscale
            tempxy[:,1]=starsFOV['Y']/fieldscale
            t3=time.time()
            #starsFOV['dmag_cloud_image']=interpolate.LinearNDInterpolator(cloudxy,np.ravel(cloudImage))(tempxy) #here's my original attempt at interpolation
            starsFOV['dmag_cloud_image']=interpolate.NearestNDInterpolator(cloudxy,np.ravel(cloudImage))(tempxy) #just use nearest neighbor.
            t4=time.time()
            #print '%f to generate clouds, %f to scale, %f to interpolate %i stars'%(t2-t1,t3-t2,t4-t3, np.size(starsFOV['X']))
            #print np.mean(cloudImage),np.std(cloudImage)
        else:
            starsFOV['dmag_cloud_image']=starsFOV['dmag_color_dist']*0

        #use Kepler stellar variability distribution
        if variability:
            star_match = np.in1d(sub_variability['id'], starsFOV['id'])
            starsFOV['dmag_kep']=variOffset(sub_variability['R_var'][star_match],sub_variability['type'][star_match])
        else:
            starsFOV['dmag_kep'] = starsFOV['rmagtrue']*0.

        #generate variable stars
        starsFOV['dmag_rr']=starsFOV['dmag_zp']*0
        if rr_fraction != 0:
            rr_stars = np.where(starsFOV['id']-np.floor(starsFOV['id']*rr_fraction)/rr_fraction == 0)
            if np.size(rr_stars) > 0:
                starsFOV['dmag_rr'][rr_stars] = rr_amplitude*(np.random.rand(np.size(rr_stars))-.5)
        #make variations based on the gain
        starsFOV['dmag_gainvar']=starsFOV['rmagtrue']*0.
        if gainvar !=0:
            #15 chips x 15 chips.  so 30 gain patches by 240 
            nx=30
            ny=240
            #x,y should run from -1*halfX to halfX
            halfX, halfY = gnomonic_project_toxy((raCen+rad_fov_deg/np.cos(decCen*deg2rad))*deg2rad,
                                             (decCen+rad_fov_deg)*deg2rad, raCen*deg2rad, decCen*deg2rad)
            ampXlen = halfX*2/nx
            ampYlen = halfY*2/ny
            ampx = np.floor((starsFOV['X']+halfX) / ampXlen)
            ampy = np.floor((starsFOV['Y']+halfY) / ampYlen)
            starsFOV['amp'] = ampy*nx+ampx
            amps = np.unique(starsFOV['amp'])
            for amp in amps:
                inamp=np.where(starsFOV['amp'] == amp)
                starsFOV['dmag_gainvar'][inamp] = starsFOV['dmag_gainvar'][inamp]-2.5*np.log10(1.+blockvisits['ampvar'][visitN][inamp]*gainvar)#-2.5*n.log10(1.+random.randn(1)*gainvar)
 
        starsFOV['dmag_exptvar'] = starsFOV['rmagtrue']*0.
        if exptvar != 0:
            starsFOV['dmag_exptvar'] = starsFOV['dmag_exptvar'] - 2.5*np.log10(1.+exptvar*np.random.randn(np.size(starsFOV['dmag_exptvar'])))
        #Use an rmag that has been altered by the patch zeropoint and clouds to calc snr-based noise, rather than the true mag.
        rmag_temp = starsFOV['rmagtrue'] + starsFOV['dmag_zp'] + starsFOV['dmag_zp_dist'] + \
                    starsFOV['dmag_cloud'] +starsFOV['dmag_cloud_image']+starsFOV['dmag_rr'] 
        if errorDist == "Cauchy":
            stats.cauchy(std=1.)
            starsFOV['dmag_var'] = stats.cauchy.rvs(size=len(starsFOV['ra']))*mag_rand_err
            starsFOV['dmag_var'] = starsFOV['dmag_var'].clip(min=-0.1,max=0.1)
            #starsFOV['dmag_snr'] = stats.cauchy.rvs(size=len(starsFOV['ra'])) * calcMagErrors(starsFOV['rmagtrue'], m5=m5, filter=filter)
            #starsFOV['dmag_snr'] = stats.cauchy.rvs(size=len(starsFOV['ra'])) * calcMagErrors(rmag_temp, m5=m5, filter=filter, error_sys=0)
            #starsFOV['dmag_snr'] = starsFOV['dmag_snr'].clip(min=-0.1,max=0.1)
        if errorDist == "Gaussian":
            # First calculate random error in this mag value - jitter or variability
            starsFOV['dmag_var'] = mag_rand_error * random.randn(len(starsFOV['ra']))
            # also calculate a magnitude based error (related to SNR)
            #starsFOV['dmag_snr'] = random.randn(len(starsFOV['ra'])) * \
            #                       calcMagErrors(starsFOV['rmagtrue'], m5=m5, filter=filter, error_sys=0.0)
        starsFOV['dmag_snr'] = random.randn(len(starsFOV['ra'])) * \
                               calcMagErrors(rmag_temp, m5=m5, filter=filter, error_sys=0.) #make snr error always Gaussian
        # so this random error includes random jitter plus SNR error
        starsFOV['dmag_rand'] = np.sqrt(starsFOV['dmag_snr']**2 + starsFOV['dmag_var']**2)
        starsFOV['rmagobs'] = (starsFOV['rmagtrue'] + starsFOV['dmag_var'] + starsFOV['dmag_snr'] + \
                              starsFOV['dmag_color'] + starsFOV['dmag_color_dist'] + \
                              starsFOV['dmag_zp'] + starsFOV['dmag_zp_dist'] + \
                              starsFOV['dmag_sin'] + starsFOV['dmag_flat_sin'] + starsFOV['dmag_cloud'] + \
                              starsFOV['dmag_cloud_image']+starsFOV['dmag_rr']+starsFOV['dmag_kep']+starsFOV['dmag_gainvar']+
                              starsFOV['dmag_exptvar']+starsFOV['dmag_illumcorr'])

        # add temperature effects on magnitude, if these have been calculated
        if starsFOV.has_key('dmag_temp'):
            starsFOV['rmagobs'] += starsFOV['dmag_temp']
            
        # Calculate reported error
        if calcerror:
            starsFOV['magerr'] = calcMagErrors(starsFOV['rmagtrue'], m5=m5, filter=filter)
        else:
            cut_off=0.001
            #temp =  starsFOV['rmagtrue']*0+np.std(starsFOV['rmagtrue']-starsFOV['rmagobs'])
            temp = starsFOV['rmagtrue']*0+np.abs(starsFOV['rmagtrue']-starsFOV['rmagobs']+starsFOV['dmag_zp'])
            #temp[np.where(temp == 0)] = np.std(starsFOV['rmagtrue']-starsFOV['rmagobs'])#mag_rand_error
            #temp[np.where(temp == 0)] = mag_rand_error
            #print np.min(temp), np.max(temp), len(np.ravel(np.where(temp == 0)))
            if len(temp) != 0:
                if np.min(temp) < cut_off:
                    temp[np.where(temp < cut_off)] = np.std(starsFOV['rmagtrue']-starsFOV['rmagobs'])
                    temp[np.where(temp < cut_off)] =  mag_rand_error
            starsFOV['magerr'] = temp#starsFOV['rmagtrue']*0+ temp#temp #starsFOV['rmagtrue']*0+ temp
            starsFOV['magerr']=starsFOV['magerr'].clip(min=systematicErr)
        
        # npatches gives you the number of patches in the entire FOV.
        # Assume they are squarely distributed.
        # Calculate patch identification; i=0,1,2,3; j=0,1,2,3.
        patch_per_side = int(np.sqrt(npatch))
        halfX, halfY = gnomonic_project_toxy((raCen+rad_fov_deg/np.cos(decCen*deg2rad))*deg2rad,
                                             (decCen+rad_fov_deg)*deg2rad, raCen*deg2rad, decCen*deg2rad)
        patchX = halfX*2/patch_per_side
        patchY = halfY*2/patch_per_side
        iPatch = np.empty(len(starsFOV['X']), dtype='int')
        jPatch = np.empty(len(starsFOV['Y']), dtype='int')
        iPatch = np.floor((starsFOV['X']+halfX) / patchX) 
        jPatch = np.floor((starsFOV['Y']+halfY) / patchY) 
        # Add a little kludge here that accounts for slight rounding errors.
        iPatch = np.where(iPatch>(patch_per_side-1), (patch_per_side-1), iPatch)
        iPatch = np.where(iPatch<0, 0, iPatch)
        jPatch = np.where(jPatch>(patch_per_side-1), (patch_per_side-1), jPatch)
        jPatch = np.where(jPatch<0, 0, jPatch)
        patchMin = 10
        starsFOV['subpatch'] = patch_per_side*iPatch + jPatch
        starsFOV['fullpatch'] = npatch*blockvisits['id'][visitN] + starsFOV['subpatch']+patchMin
        #crop off corner patches
        #badpatches=[0,patch_per_side-1, npatch-patch_per_side, npatch-1]
        #condition = np.where( (np.array(starsFOV['subpatch']) != badpatches[0]) &
        #                     (np.array(starsFOV['subpatch']) != badpatches[1]) &
        #                     (np.array(starsFOV['subpatch']) != badpatches[2]) &
        #                     (np.array(starsFOV['subpatch']) != badpatches[3]) )
        #condition = np.where(  ((starsFOV['X'] > -1.*halfX/2  ) | (starsFOV['Y'] >  -1.*halfY/2   ) ) &
        #                      ((starsFOV['X'] > -1*halfX/2   ) | (starsFOV['Y'] <  halfY/2        ) ) &
        #                      ((starsFOV['X'] < halfX/2      ) | (starsFOV['Y'] >  -1.*halfY/2   ) ) &
        #                      ((starsFOV['X'] <  halfX/2     ) | (starsFOV['Y'] <  halfY/2      )  ) )
        condition = np.where(  ((starsFOV['X'] > -1.*halfX/5*3  ) | (starsFOV['Y'] >  -1.*halfY/5*3   ) ) &
                              ((starsFOV['X'] > -1*halfX/5*3   ) | (starsFOV['Y'] <  halfY/5*3        ) ) &
                              ((starsFOV['X'] < halfX/5*3      ) | (starsFOV['Y'] >  -1.*halfY/5*3   ) ) &
                              ((starsFOV['X'] <  halfX/5*3     ) | (starsFOV['Y'] <  halfY/5*3      )  ) &
                              (starsFOV['ccd_Name'] != '') )


        for key in starsFOV.keys():
            if key != 'night': #wtf, why didn't I make this the same length as everything else?
                starsFOV[key]=starsFOV[key][condition]
        # Calculate patch data for comparison to output of selfcal solver
        patchdat = {}
        for patch in range(0, npatch, 1):
            condition = (starsFOV['subpatch'] == patch)
            # Count the number of stars in patch.
            patchdat['count'] = len(starsFOV['ra'][condition])
            if patchdat['count'] == 0:
                continue
            # Set the visit information.
            patchdat['visitid'] = blockvisits['id'][visitN]
            patchdat['fullpatch'] = npatch * blockvisits['id'][visitN] + patch+patchMin 
            patchdat['subpatch'] = patch 
            # Calculate the median values of stars within this patch. 
            patchdat['ra'] = np.median(starsFOV['ra'][condition])
            patchdat['dec'] = np.median(starsFOV['dec'][condition])
            patchdat['rmagtrue'] = np.median(starsFOV['rmagtrue'][condition])
            patchdat['color'] = np.median(starsFOV['color'][condition])
            # Calculate median of magnitude offset of stars (obs vs. true)
            temp_dmag = starsFOV['rmagobs'][condition] - starsFOV['rmagtrue'][condition]
            patchdat['dmag'] = np.median(temp_dmag)
            patchdat['dmag_rms'] = np.std(temp_dmag)
            # Calculate the random error offset info.
            temp_mag = starsFOV['dmag_rand'][condition]
            patchdat['dmag_rand'] = np.median(temp_dmag)
            patchdat['dmag_rand_rms'] = np.std(temp_dmag)
            # Calculate the zeropoint (gray) offset info.
            temp_dmag = starsFOV['dmag_zp'][condition] + starsFOV['dmag_zp_dist'][condition]
            patchdat['dmag_zp'] = np.median(temp_dmag)
            patchdat['dmag_zp_rms'] = np.std(temp_dmag)
            # Calculate the patch color offset info.
            temp_dmag = starsFOV['dmag_color'][condition] + starsFOV['dmag_color_dist'][condition]
            patchdat['dmag_color'] = np.median(temp_mag)
            patchdat['dmag_color_rms'] = np.std(temp_dmag)
            pfile = print_patchfile(pfile, patch=patchdat)
        # Print output for selfcal solver input, on all objects through visit.
        ofile = print_obsfile(ofile, stars=starsFOV, subpatch=print_subpatch)
        # Print master output for star measurements (all info), on all of visit.
        mfile = print_masterfile(mfile, visitid=blockvisits['id'][visitN], stars=starsFOV,
                                 raCen=raCen, decCen=decCen)
        # Move on to the next visit.
    return id_next


def calcMagErrors(magnitudes, filter='r', m5=24.7, error_sys=0.004):
    """Right now this is assuming airmass of 1.2 and median sky brightness for all observations.  """
    xval=magnitudes * 0.0
    error_rand = magnitudes * 0.0
    magnitude_errors = magnitudes * 0.0
    error_data={} #data from Ivezic 2008 Table 2
    error_data['filters']=np.array(['u','g','r','i','z','y'])
#    error_data['m5']=np.array([23.9,25.0,24.7,24.0,23.3,22.1])
#    error_data['delta_m5']=np.array([0.21,0.16,0.14,0.13,0.13,0.13])
    error_data['gamma']=np.array([0.037,0.038,0.039,0.039,0.040,0.040])
    matched=np.where(error_data['filters'] == filter)
    xval = np.power(10., 0.4*(magnitudes - m5))
    error_rand = np.sqrt( (0.04-error_data['gamma'][matched])*xval+error_data['gamma'][matched]*xval*xval)
    magnitude_errors = np.sqrt(error_sys*error_sys+error_rand*error_rand)
    return magnitude_errors

#def calcMagErrors(magnitudes, error_sys=_systematicErr):
#    """ This calculates the expected LSST error in magnitudes for an array of mag measurements"""
#    xval = magnitudes * 0.0
#    error_rand = magnitudes * 0.0
#    magnitude_errors = magnitudes * 0.0
#    rgamma = 0.039       # value in r band
#    xval = np.power(10, 0.4*(magnitudes - _systemM5))
#    error_rand = np.sqrt((0.04-rgamma)*xval + rgamma*xval*xval)
#    magnitude_errors = np.sqrt(error_sys*error_sys + error_rand*error_rand)     
#    return magnitude_errors

"""
fp_map_xy is passed paired arrays of x- and y-coordinates in a tangent plane.  These are scaled to mm by 
multiplying by the telescope focal length in mm.   The names of the detector (like 'R00S21C04'), and the coordinates within
the detector are returned as arrays
"""
def fp_map_xy(xArray, yArray):
    strtype = np.dtype('S9')
    n = len(xArray)
    detName = np.zeros(n, dtype=strtype)
    detIpix = np.zeros(n, dtype=np.int)
    detJpix = np.zeros(n, dtype=np.int)
    iret = cc.new_intp()
    jret = cc.new_intp()
    
    for i in range(n):
        (ierr, name) = cc.ccs2ampbrief(xArray[i]*focalLength, yArray[i]*focalLength, iret, jret)
        if not ierr:
            detName[i] =  name
            detIpix[i] = cc.intp_value(iret)
            detJpix[i] = cc.intp_value(jret)
#            print i, xArray[i]*focalLength, yArray[i]*focalLength, detName[i], detIpix[i], detJpix[i]
        else:
            detName[i] = ''
            detIpix[i] = -1
            detJpix[i] = -1
#            print i, ierr
            
    return (detName, detIpix, detJpix)

def fp_map_temp(fpTempModel, ccdName, ccdX, ccdY, time):
    """ ccdName is like R01S12C04.   Need to split into raft (R01) and ccd (S12) parts to 
        get temp.
        
        Ignore time for now.   Later shift between two different models.
        """
    n = len(ccdName)
    ccdTemp = np.zeros(n, dtype=np.float)
    
    for i in range(n):
        if ccdName[i]:
            raft = ccdName[i][0:3]
            ccd = ccdName[i][3:6]
            ccdTemp[i] = fpTempModel.getTemp(raft, ccd, ccdX[i], ccdY[i])
#            print i, ccdName[i], raft, ccd, ccdX[i], ccdY[i], ccdTemp[i]
        else:
            ccdTemp[i] = None
        
    return ccdTemp
        
    
        
def print_visitfile(vfile, visit_output="visit.dat", visits=None):
    if vfile==None:
        vfile = open(visit_output, 'w')
        print >>vfile, "%s %s %s %s %s %s %s %s %s %s %s %s %s" %("#visitID", "raCen", "decCen",
                                                         "skyRot", "time", "zpOffset",
                                                         "zpGradient", "zpGradDir",
                                                         "colorterm", "colortermRad", "filterXoff", "filterYoff", "Night")
    if visits!=None:
        for visit in np.arange(np.size(visits['id'])):
            print >>vfile, "%d %f %f %f %f %f %f %f %f %f %f %f %f" %(visits['id'][visit],
                                                             visits['ra'][visit], visits['dec'][visit],
                                                             visits['skyrot'][visit],
                                                             visits['time'][visit],
                                                             visits['zpOff'][visit],
                                                             visits['zpGrad'][visit],
                                                             visits['zpGradDir'][visit],
                                                             visits['colorterm'][visit],
                                                             visits['colortermRad'][visit],
                                                             visits['filterXoff'][visit],
                                                             visits['filterYoff'][visit],
                                                             visits['night'][visit])
    return vfile
        
def print_starstruth(sfile, stardata_output="stardata.dat", stars=None):
    # This file contains the true stellar information
    if sfile==None: 
        sfile = open(stardata_output, 'w')  
        print >>sfile, "%s %s %s %s %s %s" % ("#StarId","#StarDBId", "ra", "dec", "magTrue", "colTrue")
    if stars != None:
        # I agree - the star['id'] and star['dbid'] are somewhat redundant.
        # However, adding a star['id'] which is guaranteed sorted / increasing makes
        # some of the later analysis of the output easier. 
        for star in range(0, len(stars['id'])):
            print >>sfile, "%d %d %f %f %f %f" % (stars['id'][star], stars['dbid'][star],
                                               stars['ra'][star], stars['dec'][star],
                                               stars['rmagtrue'][star], stars['color'][star])
    return sfile

def print_obsfile(ofile, starobs_output="star_obs.dat", stars=None, subpatch=False):
    # This is the file that goes to self-calibration solver
    if ofile==None:
        ofile = open(starobs_output, 'w')
        #if subpatch:
        print >>ofile, "%s %s %s %s %s %s %s %s %s %s" %("#PatchID", "StarID", "StarObsMag", "StarMagObsErr", "ScaledRadius", "HPid1", "HPid2","HPid3","HPid4","HPid5")
        #else:
        #    print >>ofile, "%s %s %s %s" %("#PatchID", "StarID", "StarObsMag", "StarMagObsErr")
    if stars!=None:
        for obj in range(0, len(stars['id'])):
            if subpatch:
                starRadius = np.sqrt(stars['X'][obj]**2+stars['Y'][obj]**2)/fieldscale
                print >>ofile, "%d  %d  %f  %f %f %d %d %d %d %d" %(stars['fullpatch'][obj], stars['id'][obj],
                 stars['rmagobs'][obj], stars['magerr'][obj], starRadius, stars['HPid1'][obj], stars['HPid2'][obj], stars['HPid3'][obj], stars['HPid4'][obj],stars['HPid5'][obj]) 
            else:
                print >>ofile, "%d  %d  %f  %f %d %d %d %d %d %d" %(stars['fullpatch'][obj], stars['id'][obj],
                stars['rmagobs'][obj], stars['magerr'][obj], 0, stars['HPid1'][obj], stars['HPid2'][obj], stars['HPid3'][obj], stars['HPid4'][obj],stars['HPid5'][obj])
    return ofile

def print_patchfile(pfile, patch_output="patchdata.dat", patch=None):
    # This file is used for comparing output of selfcal_solver/simple.x
    if pfile == None:
        pfile = open(patch_output, 'w')
        print >>pfile, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" \
              %("#visitID", "fullpatch", "subpatch", "nstars",
                "ra", "dec", "rmagtrue", "color", "dmag", "dmag_rms",
                "dmag_rand", "dmag_rand_rms", "dmag_zp", "dmag_zp_rms",
                "dmag_color", "dmag_color_rms" )
    else:
        print >>pfile, "%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f" \
              %(patch['visitid'], patch['fullpatch'], patch['subpatch'],
                patch['count'], patch['ra'], patch['dec'],
                patch['rmagtrue'], patch['color'],
                patch['dmag'], patch['dmag_rms'],
                patch['dmag_rand'], patch['dmag_rand_rms'],
                patch['dmag_zp'], patch['dmag_zp_rms'],
                patch['dmag_color'], patch['dmag_color_rms'])
    return pfile


def print_masterfile(mfile, master_output="master_cal.dat", visitid=None, stars=None,
                     raCen=None, decCen=None):
    # This file keeps (ALL) output info for evaluation
    if mfile==None:
        mfile = open(master_output, 'w') 
        print >>mfile, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" \
              %("#VisitID", "PatchID", "subPatchID", "StarID", "StarDBID", "StarObsMag", "StarMagErr",
                "dMag_var", "dMag_snr", "dMag_zp", "dMag_zp_dist", "dMag_color", "dMag_color_dist",
                "dMag_sin", "dMag_flat_sin", "dMag_cloud", "dMag_cloud_image","dMag_rr", "StarTrueMag", "StarTrueCol", "StarRA",
                "StarDec", "RAcen_fov", "DecCen_fov",
                "StarX", "StarY", "Night", "dm5", "dMag_kep", "dMag_gainvar", "dMag_exptvar", "dMag_illumcor", "dtemp", "dmag_temp")
    if visitid!=None:
        # Print master output for star measurements (all info), on all of visit.
        for obj in range(0, len(stars['ra']), 1):
            print >>mfile, "%d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %f %f %f %f %f %f %f" \
                  % (visitid, stars['fullpatch'][obj], stars['subpatch'][obj],
                     stars['id'][obj], stars['dbid'][obj], stars['rmagobs'][obj], stars['magerr'][obj],
                     stars['dmag_var'][obj], stars['dmag_snr'][obj],
                     stars['dmag_zp'][obj], stars['dmag_zp_dist'][obj],
                     stars['dmag_color'][obj], stars['dmag_color_dist'][obj], stars['dmag_sin'][obj],
                     stars['dmag_flat_sin'][obj], stars['dmag_cloud'][obj],stars['dmag_cloud_image'][obj],
                     stars['dmag_rr'][obj],stars['rmagtrue'][obj], stars['color'][obj], stars['ra'][obj],
                     stars['dec'][obj],
                     raCen, decCen, stars['X'][obj], stars['Y'][obj], stars['night'], stars['dm5'][obj], stars['dmag_kep'][obj],
                     stars['dmag_gainvar'][obj], stars['dmag_exptvar'][obj], stars['dmag_illumcorr'][obj],
                     stars['dtemp'][obj], stars['dmag_temp'][obj])
    return mfile




################


if __name__ == "__main__":
    # Main executable. 
    print 'starting __main__'
    import sys
    if len(sys.argv) == 1 :
        parameter_filename = "NULL"
        print "No parameter file given: will use defaults"
    else:
        parameter_filename = sys.argv[1]
        print "Using parameter file %s" %(parameter_filename)
    
    # Read input parameters from input file.
    (raMin, raMax, decMin, decMax,
     use_opsim, use_opsimdither, opsimfilter,
     use_calsim, calsimtable, dbOut, 
     radius_fov, nPatch,
     nEpoch, tstart, tstop,
     mag_rand_err, calcerror,zp_var_max, zp_grad_max, colterm_max, colterm_rad_max, color_obs_noise, color_correction_error,
     filter_jitter,gainvar,exptvar, nside, healshift, use_cloudsimage,cloud_scale,sinvarx_mag, sinvary_mag, sinvarx_scale, sinvary_scale,  sinx_phase, siny_phase, sinx_angle,
     flat_sinvarx_mag, flat_sinvary_mag, flat_sinvarx_scale, flat_sinvary_scale, flat_sinx_angle, phase_persist, 
     cloud_mag, cloud_sinvar_scale, rr_fraction, rr_amplitude,kepler_variability,dith_raOff_frac, dith_decOff_frac,
     dith_skyRot_min, dith_skyRot_max,
     magmin, magmax, nStarTot, id_start, colmin, colmax, random_seed, errorDist, systematicErr, fluxNoise, fluxFrac, 
     starobs_filename, print_subpatch, master_filename, stardata_filename, patch_filename, visit_filename, illumcorr_filename,
     fpModel1_filename, fpModel2_filename, fpTempTimescale, ymagTempColorSlope, ymagTempSlope)  \
     = readParameters(parameter_filename)

    # Generate fields of view over RA/Dec boundary, with nepoch visits per fov.
    # Add dithering and rotation angles.
    if (use_opsim):
        # Opsim generates its own rotation angle distribution, but does not add its own dither.
        visits = getVisits_db(raMin, raMax, decMin, decMax, tstart, tstop, opsimfilter, 
                              radius_fov, use_opsimdither, dith_raOff_frac, dith_decOff_frac, random_seed)
    else:
        # Generate visit data in flat random distribution. 
        visits = generateVisits(raMin, raMax, decMin, decMax, nEpoch, radius_fov,
                                dith_raOff_frac, dith_decOff_frac, dith_skyRot_min, dith_skyRot_max,
                                random_seed)

    # Generate zeropoint errors, gradients, colorterms for each visit.
    visits = createVisitErrors(visits, zp_var_max, zp_grad_max, colterm_max, colterm_rad_max,
                               filter_jitter, visit_filename, random_seed, phase_persist,gainvar)
    # Calculate fieldscale - value of "X/Y" of gonomic projection at edge of field of view
    # (the point being, I would like 'X' to scale to 1 at edge, so that zp_grad_max can
    #  correspond directly to the maximum zeropoint offset expected from the gradient)
    xmax, ymax = gnomonic_project_toxy(radius_fov*deg2rad, radius_fov*deg2rad, 0, 0)
    print "# fieldscale in gnomonic projection %f %f" %(xmax, ymax)
    fieldscale = max(xmax, ymax)

    # If we're running for the 'y' filter, prepare for using temperature dependent magnitudes
        
    if opsimfilter == 'y' and fpModel1_filename:
        fpTempModel1 = fpm.FocalplaneThermalModel(fpModel1_filename)
        fpTempModel2 = fpm.FocalplaneThermalModel(fpModel2_filename)
    else:
        fpTempModel1 = None
        fpTempModel2 = None
        
    #Generate a catalog of stellar variability based on Kepler stats if requested
    if kepler_variability == True:
        variability = {'id':np.arange(nStarTot), 'R_var':np.zeros(nStarTot+0.), 'type':np.zeros(nStarTot)}
        variability['R_var']=keplerStats(nStarTot) #assign a variability magnitude
        variability['type']=variType(variability['R_var']) #assign a type of variability
    else:
        variability = False
    # Open files for output.
    sfile = print_starstruth(sfile=None, stardata_output = stardata_filename)
    ofile = print_obsfile(ofile=None, starobs_output = starobs_filename, subpatch=print_subpatch)
    pfile = print_patchfile(pfile=None, patch_output = patch_filename)
    mfile = print_masterfile(mfile=None, master_output=master_filename)
    # Connect to star database, if use_calsim.
    if (use_calsim):        
        # Connect to the database.
        conn, cursor = ui.sqlConnect(hostname=calsimhost, username=calsimusername, passwdname=calsimpassword,
                                     dbname=calsimdb, dbtype=calsimdbtype)
    else:
        cursor = None
    # Loop through sky in blocks and calculate stellar magnitudes in each block.
    #id_start = 0
    # Calculate edges of blocks.  Do these need to be even?
    rablocksize = 20 #6 #20
    decblocksize = 10 #5 #10
    if (raMin < raMax):
        rablocks = np.arange(raMin, raMax, rablocksize)
    else:
        rablocks = np.arange(raMin-360, raMax, rablocksize)
        rablocks = wrapRA(rablocks)
    decblocks = np.arange(decMin, decMax, decblocksize)
    nblocks = len(rablocks) * len(decblocks)
    if verbose:
        print "Sky divided into %d blocks." %(nblocks)
    # Calculate total sky area in simulation.
    junk = np.sin(decMax*deg2rad) - np.sin(decMin*deg2rad)
    if (raMin < raMax):
        total_sky_area = (raMax - raMin) * junk
    else:
        total_sky_area = (raMax + 360 - raMin) * junk
    for d in range(len(decblocks)):
        dec_min_block = decblocks[d]
        if d < (len(decblocks)-1):
            dec_max_block = decblocks[d+1]
        else:
            dec_max_block = decMax
        for r in range(len(rablocks)):
            ra_min_block = rablocks[r]
            if r < (len(rablocks)-1):            
                ra_max_block = rablocks[r+1]
            else:
                ra_max_block = raMax
            # Calculate how many stars per block are desired. Blocks != equal area.
            junk = np.sin(dec_max_block*deg2rad) - np.sin(dec_min_block*deg2rad)
            if ra_min_block < ra_max_block:
                block_sky_area = junk * (ra_max_block - ra_min_block)
            else:
                block_sky_area = junk * (ra_max_block + 360 - ra_min_block)
            nstars_block = np.floor(nStarTot * block_sky_area / total_sky_area)
            if verbose:
                print "Working on RA %f to %f; Dec %f to %f." %(ra_min_block, ra_max_block,
                                                               dec_min_block, dec_max_block)
                print " Generating %d stars in this block. " %(nstars_block)
                print " Area in this block is %f, total area is %f" %(block_sky_area, total_sky_area)
            id_start = generateStarMags(visits, ra_min_block, ra_max_block,
                                            dec_min_block, dec_max_block,
                                            magmin, magmax, colmin, colmax,color_obs_noise,color_correction_error,
                                            cloud_mag, cloud_sinvar_scale,rr_fraction, rr_amplitude,variability,
                                            radius_fov, fieldscale,
                                            ofile, print_subpatch, mfile, sfile, pfile, fpTempModel1, fpTempModel2,
                                            mag_rand_err, gainvar, exptvar, nside,healshift,calcerror,npatch=nPatch,
                                            use_calsim=use_calsim,use_cloudsimage=use_cloudsimage,cloud_scale=cloud_scale,
                                            cursor=cursor, calsimtable=calsimtable,
                                            nstars=nstars_block, id_start=id_start,
                                            random_seed=random_seed, fluxNoise=fluxNoise, fluxFrac=fluxFrac,
                                            filter=opsimfilter, errorDist=errorDist, systematicErr=systematicErr)

    
    # End of looping over visits, printed everything out, we're done.
    if dbOut:
        #local_conn, local_cursor = ui.sqlConnect(hostname=selfcalhost, username=selfcalusername, passwdname=selfcalpassword,
        #                                         dbname = selfcaldb, type = selfcaldbtype)
        #now make a call that just dumps the flat files into the local database
        # I think this is not implemented here - see selfcal2postgres.py instead. 
        pass
    if (use_calsim):
        ui.sqlEndConnect(conn, cursor)
    sfile.close()
    ofile.close()
    mfile.close()
    pfile.close()
