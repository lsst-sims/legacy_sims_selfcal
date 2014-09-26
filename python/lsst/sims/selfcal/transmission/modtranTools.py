"""
modtranTools.py
--
This tool box provides a set of tools needed for processing
astronomical coordinates.
Currently it contains
  GreenwichSiderealTime : to get local sidereal time at given MJD and longitude

  equatorial2local : converts target's equatorial coordinates (RA,DEC)
                     to local coordinates (azimuth, zenith angle) given
                     Universal Time (as Modified Julian Date)
                     Longitude and Latitude are fixed to LSST value
                     Times are not guaranteed to better than 0.1 second
                     between 1980 and 2030.
                     Zenith angles ignore atmospheric refraction.
                     Refer to document "Anisotropic Transmission" for details

  zenith2airmass : from z_angle to airmass : Hardie (1962) approximation
                   for sea level

  airmass2zenith : from z_angle to airmass (polynomial approximation fitted
                   to Modtran4 numerical integration at Cerro Pachon lattitude
                   and altitude)
"""

import numpy

# conversion degree to radian
D2R = numpy.pi / 180.
R2D = 1. / D2R

# January 1st 2000
Y2K_MJD = 51544.5         # in mjd
Y2K_GMST = 18.697374558   # Greenwich mean sidereal time (hours)

# LSST coordinates (decimal degrees)
LSSTLONG = -70.749389
LSSTLAT = -30.244333

# fitted parameters to compute airmass at Cerro Pachon
AM_ADD = 0.00003873
AM_NORM = 1.00003873
AM_A = 0.00548117
AM_B = 0.00832316
AM_C = 0.00711221


def GreenwichSiderealTime(mjd):
    """Return Greenwich apparent sidereal time given MJD"""
    # time (UT) ellapsed since January 1.0 2000
    delta_T = mjd - Y2K_MJD
    # Greenwich mean sidereal time (hours)
    #  (approximation according to USNO)
    Gmst = Y2K_GMST + 24.06570982441908 * delta_T
    ## Computing nutation in longitude (dpsi)
    # longitude of ascending node(deg)
    omega = 125.04 - 0.052954 * delta_T
    # mean longitude of the Sun ( " )
    L = 280.47 + 0.98565 * delta_T
    dpsi = (- 0.000319 * numpy.sin(D2R*omega) -
            0.000024 * numpy.sin(2.*D2R*L))
    ## correct Gmst to apparent by adding equation of
    ##  equinox (dpsi*cos(epsilon))
    # obliquity of equator (deg)
    epsilon = 23.4393 - 0.0000004 * delta_T
    Gst = Gmst + dpsi * numpy.cos(D2R*epsilon)
    return numpy.mod(Gst, 24.)


def equatorial2local(ra, dec, mjd, unit='deg'):
    """Convert equatorial coordinates (RA,DEC in 'unit')
    of an astronomical target to local coordinates (azimuth,
    zenith angle) given observation time (mjd), observer's
    longitude and latitude
    INPUT : ra_target, dec_target, mjd_observation
    OUTPUT : azimuth, z_angle in 'unit'
    """
    # get local sidereal time
    lst = GreenwichSiderealTime(mjd) + LSSTLONG/15.
    if unit == 'deg':
        dec = dec * D2R
    elif unit == 'rad':
        ra = ra * R2D
    else:
        raise AttributeError('Unit unknown !')
    # get hour angle (degrees)
    H = lst * 15. - ra
    Hr = H * D2R
    # trigo
    cHr, sHr = numpy.cos(Hr), numpy.sin(Hr)
    cdec, sdec = numpy.cos(dec), numpy.sin(dec)
    clat, slat = numpy.cos(LSSTLAT*D2R), numpy.sin(LSSTLAT*D2R)
    # get local coordinate projections (azimuth and zenith angle)
    cz = slat * sdec + clat * cdec * cHr
    z_angle = numpy.arccos(cz)
    #
    szsa = cdec * sHr
    szca = -1.0 * clat * sdec + slat * cdec * cHr
    szi = 1. / numpy.sqrt(1. - cz**2)
    sa = szsa * szi
    ca = szca * szi
    azimuth = numpy.arctan2(sa, ca)
    if unit == 'dec':
        return azimuth * R2D, z_angle * R2D
    elif unit == 'rad':
        return azimuth, z_angle
    else:
        raise AttributeError('Unit unknown !')


def zenith2airmass(z_angle, site='lsst', unit='deg'):
    """Compute airmass at a given zenith angle.
    INPUT options =>
     - site => 'lsst' : approximation of airmass along the slant path
                        fitted to modtran atmosphere integration above
                        the Cerro Pachon site
            => 'sea'  : compute airmass using the general purpose approximation
                        by Hardie (1962)
     - unit => 'deg' if input is in degrees
            => 'rad' -------------- radians
    """
    if unit == 'deg':
        zang_rad = z_angle * D2R
    elif unit == 'rad':
        zang_rad = z_angle
    else:
        raise AttributeError('Unit unknown !')
    secz = 1. / numpy.cos(zang_rad)
    if site == 'lsst':
        chi = numpy.log(secz)
        mdtairmass = (AM_NORM * secz - AM_ADD - AM_A * chi ** 2 +
                      AM_B * chi ** 3 - AM_C * chi ** 3)
        return mdtairmass
    elif site == 'sea':
        sz1 = secz - 1.
        airmass = (secz - 0.0018167 * sz1 - 0.002875 * sz1**2 -
                   0.0008083 * sz1**3)
        return airmass


def airmass2zenith(airmass, unit='deg'):
    """Compute zenith angle given airmass at Cerro Pachon site
    mdtram2za(mdtram(z_angle)) = z_angle
    Option :
     - unit => 'deg' if output wanted in degrees
            => 'rad' ------------------- radians
    """
    sz = (airmass + AM_ADD) / AM_NORM
    chi = numpy.log(sz)
    secz = sz + AM_A * chi**2 - AM_B * chi**3 + AM_C * chi**3
    z_angle = numpy.arccos(1./secz)
    if unit == 'deg':
        return z_angle * R2D
    elif unit == 'rad':
        return z_angle
    else:
        raise AttributeError('Unit unknown !')


def movingAverage(data, window_size):
    """Compute the moving average of a series of points
    given a window size"""
    # make a loop on the data
    i = window_size
    interval = numpy.concatenate([data[-i:], data, data[:i]])
    window = numpy.ones(int(window_size))/float(window_size)
    return numpy.convolve(interval, window, 'same')[i:-i]
