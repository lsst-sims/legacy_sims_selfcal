#===========================================================================
#
#               Module : Astrotools
#
#    This tool box provides basic tools needed for processing 
#       astronomical coordinates;
#               currently (03-2010)it contains 
#
#       tsid  : to get local sidereal time at  
#               given MJD and longitude
#       eq2loc : converts target's equatorial coordinates
#                       (RA,DEC) 
#               to local coordinates
#                       (azimuth, zenith angle)
#               given Universal Time, ( given as Modified Julian Date )
#                      Longitude and Latitude are fixed to LSST value
#               Times are not guaranteed to better than 0.1 second
#                       between 1980 and 2030.
#               Zenith angles ignore atmospheric refraction.
#                       Refer to document
#                               "Anisotropic Transmission" for details  
#       airmass : from z_angle to airmass : Hardie (1962) approximation.
#
#       mdtram : from z_angle to airmass (polynomial approximation fitted
#               to Modtran4 numerical integration at Cerro Pachon lattitude
#               and altitude)
#
#---------------------------------------------------------------------------
#

import numpy
from numpy import *

#
# conversion from degree to radian
#
d2r = numpy.pi/180
r2d = 1./d2r
#
#_____________________________________________________________
#
#                       tsid
#
#    get greenwhich apparent sidereal time from given time (MJD)
#
def tsid(mjd):
        # time (UT  )ellapsed since January 1.0 2000

        delta_T = mjd - 51544.5
        
        # greenwhich mean sidereal time (hours)
        #        (approximation according to USNO)

        Tgms =  18.697374558 + 24.06570982441908 * delta_T

        # nutation in longitude (dpsi)
        omega = 125.04 - 0.052954*delta_T     #longitude of ascending node(deg)
        L = 280.47 + 0.98565*delta_T          # mean longitude of the Sun ( " )
        epsilon = 23.4393 - 0.0000004*delta_T # obliquity of equator      ( " )
        dpsi = -0.000319 * sin(d2r*omega) - 0.000024 * sin(2.*d2r*L )
        # correct Tgms to apparent by adding
        #      equation of equinox (dpsi*cos(epsilon))
        Tgs = Tgms + dpsi * cos( d2r * epsilon )
        Tgs = mod(Tgs,24.)

        # return greenwhich apparent sidereal time

        return Tgs

#________________________________________________________________
#
#                      eq2loc 
#
#     convert equatorial coordinates (RA,DEC in degrees) of an astronomical
#     target to local coordinates (azimuth, zenith angle)
#     given observation time (mjd), observer's longitude and latitude
#     
#      Input : aldeti=[alfa_target,delta_target,mjd_observation]
#              requires module tsid
#      Output : locco = [azimuth, z_angle] in decimal degrees
#
#
def eq2loc(aldeti):
        
        # fix observation coordinates (decimal degrees)
        obslong = -70.749389 # LSST longitude 
        obslat  = -30.244333 # LSST latitude
        
        # get observation data
        alfa = aldeti[0] # Right Ascencion (degrees)
        delta = aldeti[1]# Declination (degrees)
        mjd = aldeti[2]   # modified Julian date (days)
        
        # get local sidereal time
        Tls = tsid(mjd) + obslong/15.
        
        # get hour angle (degrees)
        H = (Tls*15. - alfa)

        # convert angles to radians
        Hr = H * d2r
        der = delta * d2r
        fir = obslat * d2r
        
        # get local coordinate projections ( azimuth and zenith angle)
        cz = sin(fir)*sin(der)+cos(fir)*cos(der)*cos(Hr)
        szsa = cos(der)*sin(Hr)
        szca = -cos(fir)*sin(der)+sin(fir)*cos(der)*cos(Hr)
        
        # separate trigo fun
        szi = 1./sqrt(1.-cz*cz)
        sa = szsa * szi
        ca = szca * szi
        
        # extract angles
        z_angle = arccos(cz)*r2d
        azimuth = arctan2(sa,ca)* r2d
        
        locco = [azimuth,z_angle]

        return locco
#_____________________________________________________________________________
#
#                       airmass
#
#  Compute airmass at given zenith angle (degrees)
#   using the general purpose approximation by Hardie (1962) 
#       sz1 = secz - 1
#       airmass = secz - 0.0018167*sz1 - 0.002875*sz1^2 - 0.0008083*sz1^3
#       Warning : This model is for sea level
#
#
def airmass(z_angle) :
        secz = 1./cos(z_angle*d2r)
        sz1 = secz-1.
        airmass = secz - 0.0018167*sz1 - 0.002875*sz1**2 - 0.0008083*sz1**3

        return airmass
#_____________________________________________________________________________
#
#                       mdtram (modtran airmass at Cerro Pachon)
#
# Compute airmass at given zenith angle (degrees).
#   Approximation of airmass along the slant path
#   fitted to modtran atmosphere integration above the LSST Cerro Pachon site
#   (This includes latitude and altitude)
#   z_angle must be in decimal degrees
#
#
def mdtram(z_angle):
        secz = 1./cos(z_angle*d2r)
        chi = log(secz)  
        mdtairmass =1.00003873*secz-0.00548117*chi**2\
            + 0.00832316*chi**3 - 0.00711221*chi**3 - 0.00003873
        return mdtairmass

def mdtram2za(airmass):
        sz = (airmass + 0.00003873)/1.00003873
        chi = numpy.log(sz)
        secz = sz + 0.00548117*chi**2 - 0.00832316*chi**3 + 0.00711221*chi**3
        chi = numpy.log(secz)
        secz = sz + 0.00548117*chi**2 - 0.00832316*chi**3 + 0.00711221*chi**3
        z_angle = r2d*numpy.arccos(1./secz)

        return z_angle
