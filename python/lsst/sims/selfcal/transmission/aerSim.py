#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import numpy as np
# import matplotlib.pyplot as plt
import numpy.fft as npf
import numpy.random as npr
import numpy.linalg as npl
from scipy.interpolate import InterpolatedUnivariateSpline  # , splev, splrep

# LOCAL
#------
# mainpath = '/Users/alexandreboucaud/work/LSST/calib/aerosols/'
# plotpath = mainpath + 'plots/'

# LSST
#-----
mainpath = os.getenv('ATMOSPHERE_TRANSMISSION_DIR')

mauna = os.path.join(mainpath, 'data/MaunaLoa_1997-2009_final.dat')
casleo = os.path.join(mainpath, 'data/Casleo_2011-2012_final.dat')

# days in a year
YLEN = 365
# years of data
N_YEARS = 13
# number of days for output
N_DAYS = N_YEARS * YLEN
# size of the running average (in days)
N_RUN = 60
# polynome degree for fit
P_DEG = 2

# Defining absissae for plot accross full LSST spec range.
# modtran mesh is 5212 iso-energy bins
#   starting at frequency 9040 cm-1 (1106.19469027nm)
#   ending at frequency 35095 cm-1 (284.940874768nm)
NI = 500
NU_MIN = 9040.
NU_MAX = 35100
NU_STEP = (NU_MAX - NU_MIN) / NI

##########
# METHODS
##########


# GLOBAL
#-------

def scale(data):
    """Empiric relation to scale data between Hawaii and
    South America latitude"""
    return data * 1.45 + 0.004


def dictCut(dictionary, cut):
    """Apply a cut to all the dictionary items"""
    for name, arr in dictionary.items():
        dictionary[name] = arr[cut]


def movingAverage(interval, window_size):
    """Smooth an array using a convolution"""
    edge = window_size
    interval = np.concatenate([interval[-edge:], interval, interval[:edge]])
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(interval, window, 'same')[edge:-edge]


# DICTIONARY CONSTRUCTION
#------------------------

def select_daily_values(wl, ambin):
    """
    Compute the average daily value, if any.

    @param str    wl | wavelength
    @param str ambin | airmass bin number

    return dict
        array int     mjd_avg | days with at least a value
        array float  laer_avg | existing daily values
        array float laer_days | indexed daily values
        array int    mis_vals | days with no value
    """
    data_days = np.zeros(N_DAYS)
    mjd_int = np.array(mjd_dict[wl][ambin], dtype='int')
    mjd_avg = []
    data_avg = []
    mis_vals = []
    for d in xrange(N_DAYS):
        daily_values = np.where(mjd_int == d)[0]
        if daily_values.shape[0]:
            data = laer_dict[wl][ambin]
            daily_data = data[daily_values]
            avg_data = daily_data.mean()
            mjd_avg.append(d)
            data_avg.append(avg_data)
            # draw random value
            npr.shuffle(daily_data)
            data_days[d] = daily_data[0]
        else:
            mis_vals.append(d)
    tdict = {}
    tdict['mjd_avg'] = np.array(mjd_avg, dtype='int')
    tdict['laer_avg'] = np.array(data_avg, dtype='float')
    tdict['laer_days'] = data_days
    tdict['mis_vals'] = np.array(mis_vals, dtype='int')
    return tdict


def average_n_variance(tdict):
    """
    Compute the average daily value and variance
    over the 13 years of data

    @param dict tdict | dictionary issued from select_daily_values

    return dict
        array float laer_days | indexed daily values
        array int    mis_vals | days with no value
        array float laer_yavg | values averaged o/ 13y
        array float laer_ystd | std for each day o/ 13y
    """
    idx = tdict['mjd_avg'] % YLEN
    data_yavg = np.zeros(YLEN)
    data_ystd = np.zeros(YLEN)
    for d in xrange(YLEN):
        dvals = tdict['laer_avg'][np.where(idx == d)]
        davg = np.mean(dvals)
        data_yavg[d] = davg
        data_ystd[d] = np.sqrt(np.mean((dvals - davg)**2))
    nul = np.where(data_ystd == 0.0)[0]
    data_ystd[nul] = 0.5 * (data_ystd[nul + 1] + data_ystd[nul - 1])
    tdict['laer_yavg'] = data_yavg
    tdict['laer_ystd'] = data_ystd
    del(tdict['mjd_avg'])
    del(tdict['laer_avg'])
    return tdict


def fill_blanks(tdict):
    """
    Fill the empty daily values using average and variance

    @param dict tdict | dictionary issued from average_n_variance

    return dict
        array float laer_days | indexed daily values
        array float laer_yavg | values averaged o/ 13y
    """
    for idx in tdict['mis_vals']:
        y_id = idx % YLEN
        mu = tdict['laer_yavg'][y_id]
        sigma = tdict['laer_ystd'][y_id]
        tdict['laer_days'][idx] = np.random.normal(mu, sigma / 2.)
    del(tdict['mis_vals'])
    del(tdict['laer_ystd'])
    return tdict


def correct_seasonal_var(tdict, size):
    """
    Moving average and correction of data

    @param dict tdict | dictionary issued from fill_blanks
    @param int   size | window size for the running average

    return dict
        array float   laer_days | indexed daily values
        array float  laer_final | daily values corrected from
                                  seasonal variation
        interp spline laer_yspl | interpolated seasonal variation
    """
    days_arr = np.arange(N_DAYS)
    tdict['laer_yspl'] = InterpolatedUnivariateSpline(
        np.arange(YLEN), movingAverage(tdict['laer_yavg'], size))
    tdict['laer_final'] = (
        tdict['laer_days'][days_arr] -
        tdict['laer_yspl'](days_arr % YLEN))
    del(tdict['laer_yavg'])
    return tdict


def get_dict(wavelength, airmass_bin, size=N_RUN):
    """Sums up a few methods in a row, cf. to their docstrings for info"""
    tdict = select_daily_values(wavelength, airmass_bin)
    tdict = average_n_variance(tdict)
    tdict = fill_blanks(tdict)
    fdict = correct_seasonal_var(tdict, size)
    return fdict


# FFT AND RANDOMIZATION
#----------------------

def vectorize(aerlist, amlist, name='laer_final'):
    """
    Build a vector for all the combinations of a given parameter

    @param str      name | name of the wanted dictionary item

    return array
    """
    vector = []
    for wl in aerlist:
        for iam in amlist:
            data = dailydict[wl][iam][name]
            fftdata = npf.rfft(data)
            nfreq = len(fftdata) - 1
            # avgfft = fftdata.mean()
            # print wl, iam, 'Average value fft\t', avgfft
            vector.append(fftdata)
            # vector.append(avgfft)
    vector = np.asarray(vector, dtype='complex')
    return vector, nfreq


def vector2dict(vector, log=True):
    """
    Take the output vector of the randomization and work it to return
    a dictionary of readable and directly usable content.

    @param array    vector | output of randomization
    @keyword bool      log | set to True if outdata in log

    return dict
        str arr      aerstr | wavelengths of optical depth
            str arr   amstr | bin number of airmass
                array float | simulated values of aerosol optical depth data
                              with restored seasonal var (no more in log !!!)
    """
    ndict = {}
    ipar = 0
    for wl in aerstr:
        ndict[wl] = {}
        for am in amstr:
            simdata = vector[ipar, :]
            # restore the seasonal variations
            days = np.arange(len(simdata))
            spline = dailydict[wl][am]['laer_yspl']
            logdata = simdata + spline(days % YLEN)
            if log:
                ndict[wl][am] = logdata
            else:
                ndict[wl][am] = np.exp(logdata)
            ipar += 1
    return ndict


def covariance(vector):
    """"Compute the covariance matrix between optical depth
    at the five wavelength"""
    npar, nfreq = vector.shape
    avgvect = vector.mean(axis=1)
    sqavg = np.outer(avgvect, np.conj(avgvect))
    bigmatrix = np.zeros((npar, npar, nfreq), dtype='complex')
    for freq in xrange(nfreq):
        bigmatrix[:, :, freq] = np.outer(vector[:, freq], np.conj(vector[:, freq]))
    avgsq = bigmatrix.mean(axis=2)
    covmatrix = avgsq - sqavg
    return covmatrix


def randomize(eigenvals, eigenvect, seed=10):
    """Get eigenvalues and eigenvectors of the covariance matrix
    and draw the random values

    FYI: D = P-1 A P
    """
    evlen = eigenvals.shape
    npr.seed(seed)
    s = npr.randint(0, 10000, evlen)
    randvalues = np.zeros(evlen)
    #randvector = np.zeros(eigenvals.shape, dtype='complex')
    for ival, val in enumerate(eigenvals):
        npr.seed(s[ival])
        randvalues[ival] = npr.normal(0, np.sqrt(np.abs(val)))
        # npr.seed(s)
        # phi = npr.uniform(0, 2 * np.pi)
        # randvalues[ival] = amp * (np.cos(phi) + 1j * np.sin(phi))
    randvector = np.dot(eigenvect, randvalues)
    return randvector


# FITTING
#--------

def main(seed=10):
    master_vect, nfreq = vectorize(aerstr, amstr)
    # print master_vect
    npar = master_vect.shape[0]
    #print 'Vector shape: %d x %d' % (npar, nfreq)
    master_rand = np.zeros((npar, nfreq + 1), dtype='complex')
    master_out = np.zeros((npar, nfreq * 2))
    covariance_matrix = covariance(master_vect)
    # covariance_matrix = np.outer(master_vect, np.conj(master_vect))
    eigenvalues, eigenvectors = npl.eigh(covariance_matrix)
    # print 'Covariance Matrix:\n', covariance_matrix
    # print 'Eigenvalues:\n', eigenvalues
    # print 'Eigenvectors:\n', eigenvectors
    npr.seed(seed)
    sprseed = npr.randint(0, 100000, nfreq)
    for freq in xrange(nfreq):
        master_rand[:, freq] = randomize(eigenvalues, eigenvectors, sprseed[freq])
    for par in xrange(npar):
        master_out[par] = npf.irfft(master_rand[par])
    out_dict = vector2dict(master_out, log=True)
    return out_dict


def getAerParameters(seed=10, airmass='0'):
    """Return a vector containing the daily fitting parameters for aerosols
    at Cerro Pachon for 13 years

    bonus: write a temporary file containing the data"""
    if type(airmass) == int:
        airmass = str(airmass)
    # Get simulated parameters in a dictionary
    outdict = main(seed)
    # Cast them into a vector
    outvect = np.zeros((N_DAYS - 1, len(aerstr)))
    for iwl, wl in enumerate(aerstr):
        outvect[:, iwl] = outdict[wl][airmass]
    # Construct the Vandermonde matrix for standard deviation on the datapoints
    vx = np.vander(laer_wl, P_DEG + 1)
    # Create temporary file to store the results
    tmpfile = open('tempfile.dat', 'w')
    strhead = '#Time [day]\tAirmass bin\tp0\tp1\tp2\tstdev\n\
        #Seed used for this simulation = {0:d}\n'.format(seed)
    tmpfile.write(strhead)
    # initialize vector for fitting parameters
    fitvect = np.zeros((N_DAYS - 1, 4))
    for day in xrange(N_DAYS - 1):
        lvect = outvect[day]
        # Fit
        pfit = np.polyfit(laer_wl, lvect, P_DEG)
        # Standard deviation
        stdev = np.std(lvect - np.dot(vx, pfit))
        # Store in vector
        fitvect[day] = np.array([pfit[0], pfit[1], pfit[2], stdev])
        # Write in file
        strline = '{0:d}\t{1:s}\t{2:.6f}\t{3:.6f}\t{4:.6f}\t{0:.6f}\n'.format(
            day, airmass, pfit[0], pfit[1], pfit[2], stdev)
        tmpfile.write(strline)
    tmpfile.close()
    return fitvect


# PLOTTING
#---------

# def plot_spectrum(wl, ambin):
#     data = npf.rfft(dailydict[str(wl)][str(ambin)]['laer_final'])
#     k = npf.fftfreq(N_DAYS)
#     xplot = npf.fftshift(1. / k)
#     xplot = xplot[2372:]
#     fnorm = data / N_DAYS
#     fplot = np.abs(npf.fftshift(fnorm))
#     print len(xplot)
#     print len(fplot)
#     plt.plot(xplot, fplot)
#     plt.xlim(0, 200)
#     plt.ylim(0, 0.04)
#     plt.show()


# def plot_histograms(fdict, pltbins):
#     f, ax = plt.subplots(3, 2, sharex=True, sharey=True)
#     for idx, wl in enumerate(aerstr):
#         ix, iy = idx % 3, idx / 3
#         dat = np.array([])
#         sim = np.array([])
#         for am in amstr:
#             dat = np.hstack((dat, np.exp(dailydict[wl][am]['laer_days'])))
#             sim = np.hstack((sim, fdict[wl][am]))
#         ax[ix, iy].hist(dat, normed=True, bins=pltbins, alpha=0.4,
#                         histtype='stepfilled', label='original data')
#         ax[ix, iy].hist(sim, normed=True, bins=pltbins, alpha=0.4,
#                         histtype='stepfilled', label='simulated data')
#         ax[ix, iy].set_title(r'$\tau_{%s}$' % wl)
#     ax[1, 1].legend(loc=1)
#     # plt.legend(('original data', 'simulated data'), loc=1)
#     plt.xlim(0, 0.1)
#     plt.show()


# def plot_data_vs_daily(pltbins):
#     f, ax = plt.subplots(3, 2, sharex=True, sharey=True)
#     for idx, wl in enumerate(aerstr):
#         ix, iy = idx % 3, idx / 3
#         dat = np.array([])
#         for am in amstr:
#             dat = np.hstack((dat, np.exp(dailydict[wl][am]['laer_days'])))
#         ax[ix, iy].hist(aerlist[idx], normed=True, bins=pltbins, alpha=0.4,
#                         histtype='stepfilled', label='full data')
#         ax[ix, iy].hist(dat, normed=True, bins=pltbins, alpha=0.4,
#                         histtype='stepfilled', label='daily data')
#         ax[ix, iy].set_title(r'$\tau_{%s}$' % wl)
#     ax[1, 1].legend(loc=1)
#     # plt.legend(('original data', 'simulated data'), loc=1)
#     plt.xlim(0, 0.2)
#     plt.show()


# I/O
#----

# def write_daily_values():
#     amlist = ['2.0', '3.0', '4.0']
#     headerlist = ['time', 'airmass', 'tau_870', 'tau_675',
#                   'tau_500', 'tau_440', 'tau_380']
#     header = '\t'.join(h for h in headerlist)
#     with open('daily_aerosols.dat', 'w') as daer:
#         daer.write('%s\n' % header)
#         for day in xrange(N_DAYS):
#             for am in xrange(3):
#                 line = '\t'.join('%.4f' % np.exp(
#                     dailydict[wl][str(am)]['laer_days'][day])
#                     for wl in aerstr[::-1])
#                 daer.write('%d\t%s\t%s\n' % (day, amlist[am], line))
#     print 'DONE'


# def write_simulated_values(simul, finaldict):
#     amlist = ['2.0', '3.0', '4.0']
#     headerlist = ['time', 'airmass', 'tau_870', 'tau_675',
#                   'tau_500', 'tau_440', 'tau_380']
#     header = '\t'.join(h for h in headerlist)
#     with open('simulated_aerosols.dat', 'w') as faer:
#         faer.write('%s\n' % header)
#         for day in xrange(simul.shape[1]):
#             for am in xrange(3):
#                 line = '\t'.join('%.4f' % finaldict[wl][str(am)][day]
#                                  for wl in aerstr[::-1])
#                 faer.write('%d\t%s\t%s\n' % (day, amlist[am], line))
#     print 'DONE'


# def write_scaled_data():
#     headerlist = ['time', 'airmass', 'tau_870', 'tau_675',
#                   'tau_500', 'tau_440', 'tau_380']
#     header = '\t'.join(h for h in headerlist)
#     with open('casleolike_aerosols.dat', 'w') as faer:
#         faer.write('%s\n' % header)
#         for day in xrange(N_MJD):
#             line = '\t'.join('%.4f' % np.exp(laer_masterdict[wl][day])
#                              for wl in aerstr[::-1])
#             faer.write('%.3f\t%.2f\t%s\n' % (maindict['mjd'][day], maindict['airmass'][day], line))
#     print 'DONE'


####################
# BEGINNING OF FILE
####################

# Mauna Loa data
data_mauna = np.loadtxt(mauna)
mjds, mjdy, airmass, aer380, aer440, aer500, aer675, aer870, ang_e, fmf = data_mauna.T
# starts at 0
mjd = mjds - mjds[0]
# cut redondant values
dmjd = np.diff(mjd)
dmjd_cut = np.where(dmjd != 0.0)

# Main parameters dictionary
mainstr = ['mjd', 'mjdy', 'airmass', 'ang', 'fmf']
mainlist = [mjd, mjdy, airmass, ang_e, fmf]
maindict = dict(zip(mainstr, mainlist))

# Optical depth wavelengths
aerstr = ['380', '440', '500', '675', '870']
# Optical depth dictionary
aerlist = [aer380, aer440, aer500, aer675, aer870]
# rescale from Hawaii to Argentina
saerlist = [scale(aer) for aer in aerlist]
# convert to log
laerlist = [np.log(aer) for aer in saerlist]
# put in a dictionary
laer_masterdict = dict(zip(aerstr, laerlist))
# wavelength in log
laer_wl = np.array([np.log(int(aer)) for aer in aerstr], dtype='float')

# for i, arr in enumerate(aerlist):
    # print 'Wavelength:\t{0}'.format(aerstr[i])
    # print 'Length:\t{0}'.format(arr.shape[0])
    # print 'Missing values:\t{0}\n'.format(arr[np.where(arr < 0.0)].shape)

# Apply cut on dictionaries
dictCut(maindict, dmjd_cut)
dictCut(laer_masterdict, dmjd_cut)

# Number of values after cut
N_MJD = len(maindict['mjd'])

# Airmass cuts
# empiric bins created with equal data in each
am_bins = [1.0, 2.6, 3.52, 5.0]
am_cuts = []
for imin, imax in zip(am_bins[:-1], am_bins[1:]):
    am_cuts.append(np.where((maindict['airmass'] >= imin) &
                            (maindict['airmass'] < imax)))
amstr = [str(i) for i in range(len(am_cuts))]

# Dictionary with wavelengths AND airmass cuts
laer_dict = {}
mjd_dict = {}
for wl in aerstr:
    laer_dict[wl] = {}
    mjd_dict[wl] = {}
    for iam, cut in zip(amstr, am_cuts):
        laer_dict[wl][iam] = laer_masterdict[wl][cut]
        mjd_dict[wl][iam] = maindict['mjd'][cut]

# Create daily dictionaries
dailydict = {}
for wl in aerstr:
    dailydict[wl] = {}
    for iam in amstr:
        dailydict[wl][iam] = get_dict(wl, iam)
# dictlist = ['laer_days', 'laer_yspl', 'laer_final']

###########
# TESTING #
###########

# nval = 30
# imin = 650
# imax = imin + nval
# outdict = main()
# am = '0'
# numesh = np.arange(NU_MIN, NU_MAX, NU_STEP)
# lsstwl = 1.e7 / numesh  # convert frequencies (cm-1)into wvlength (nm)
# lwv = np.log(lsstwl)

# vcur = np.vander(lwv, P_DEG + 1)

# testarr = np.zeros((nval, len(aerstr)))
# for iwl, wl in enumerate(aerstr):
#     # testarr[:, iwl] = outdict[wl][am][imin:imax]
#     testarr[:, iwl] = outdict[wl][am]

# Construct the Vandermonde matrix for standard deviation
# vx = np.vander(laer_wl, P_DEG + 1)
# aer_wl = np.exp(laer_wl)

# col = ['b', 'r', 'm', 'g', 'c']
# xplot = np.arange(N_DAYS - 1)
# xplot2 = np.arange(0, N_DAYS - 1, 0.5)

# f, ax = plt.subplots(3, 2, sharex=True)
# pf = np.zeros((nval, 3))
# for day in xrange(nval):
#     lvect = testarr[day]
#     # Fit
#     pfit = np.polyfit(laer_wl, lvect, P_DEG)
#     pf[day] = pfit

# p0 = pf[:, 0]
# p1 = pf[:, 1]
# p2 = pf[:, 2]
# p0sp = splrep(xplot, p0)
# p1sp = splrep(xplot, p1)
# p2sp = splrep(xplot, p2)
# p02 = splev(xplot2, p0sp)
# p12 = splev(xplot2, p1sp)
# p22 = splev(xplot2, p2sp)

# ax[0, 0].plot(p0, 'r,')
# ax[1, 0].plot(p1, 'b,')
# ax[2, 0].plot(p2, 'g,')
# ax[0, 1].plot((p02-p0)/p0, 'r,')
# ax[1, 1].plot((p12-p1)/p1, 'b,')
# ax[2, 1].plot((p22-p2)/p2, 'g,')

# plt.show()
