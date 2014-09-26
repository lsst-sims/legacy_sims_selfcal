# -*- coding: cp1252 -*-
#==============================================================================
#
#               module makesequence.py  (alias : Pointing.py):
#
#-------------------------------------------------------------------------
# Read description of a sequence of LSST pointings  from "sequencename"_Descript.txt
# Read from Obs history a sequence of LSST Pointing according to description
#        with parameters as described in description file
# Search for relevant Atmosphere data in Snn.nnAtm_history.dat
#_______________________________________________________________________________
# common items

#       modules
from numpy import *
import Astrotools
import copy

#          constants
d2r = pi/180.
r2d = 1./d2r

#          paths
#                    Give below absolute path to mkatmos in principle it is the only change
localpath =  '/Users/ljones/lsst_home/Darwin/sims/atmosphere/'

mlsstpath = localpath +'mkatmos/SpData/'
runpath = mlsstpath + 'Runparm/'

meteosuffix = '.02'
obshistname = 'Opsim3.61.sel'
expname = 'Opsim3.61'+meteosuffix+'sel'

LSSTdatafile = 'LsstObsHist/' +obshistname+'.dat'
Atfile = 'AtmosHist/meteo'+meteosuffix+'Atm_history.dat'
Atmos_parmfile= mlsstpath + Atfile
#____________
#            Read Atmosphere history for period in list aparmlis
#



aparmf = open(Atmos_parmfile,'r')
aparmlis = list(aparmf)
aparmf.close()
#
#               
#               Adding atmosphere parameters at appropriate MJD time 
#
print 'succesfully read ', len(aparmlis),'atmos parms in file : ',Atfile

#
# extract MJD start 
#
MJDS = float( aparmlis[0].split('$')[0].split('=')[1])
mjdstep = 0.01
#
# extract Aerosol drivers
#
natmos = len(aparmlis)
avpl = [[]]* natmos # list of aerosol parameters
atpl = [[]]* natmos # list of strings containing other atmosphere parameters
for iat in range(0,natmos):
    for parv in aparmlis[iat].split('$') :
        if parv.split('=')[0].strip()=='VIS0':
            vis0 = float(parv.split('=')[1])
        if parv.split('=')[0].strip()=='VISAMP':
            visamp = float(parv.split('=')[1])
        if parv.split('=')[0].strip()=='VISAZ':
            visaz = float(parv.split('=')[1])
        if parv.split('=')[0].strip()=='ZAER11':
            zaer11 = float(parv.split('=')[1])
        if parv.split('=')[0].strip()=='ZAER12':
            zaer12 = float(parv.split('=')[1])
        if parv.split('=')[0].strip()=='SCALE1':
            scale1 = float(parv.split('=')[1])
        if parv.split('=')[0].strip()=='ZAER21':
            zaer21 = float(parv.split('=')[1])
        if parv.split('=')[0].strip()=='ZAER22':
            zaer22 = float(parv.split('=')[1])
        if parv.split('=')[0].strip()=='SCALE2':
            scale2 = float(parv.split('=')[1])
    avpl[iat] = [vis0,visamp,visaz,zaer11,zaer12,scale1,zaer21,zaer22,scale2]
    atpl[iat] = [aparmlis[iat].split('VIS0')[0]]
print 'Aerosol drivers set for ' , iat, 'atmospheres'
print 'last are :' , avpl[iat],atpl[iat]
#

#_______________________________________________________________________________
#
#               Observing History acquisition
#               including MJD start and MJD end
#

mydata = mlsstpath + LSSTdatafile
#
# open obs history file and convert to list
#
mydatf = open(mydata,'r')
fullhlines= list(mydatf)
mydatf.close()

#
# get list of retrievable parameters in file
headers = fullhlines[0]
head = headers.split('\t')
for i in range(0,len(head)):
	head[i]=head[i].strip()
#  list of obs.history parameters to be extracted
headextra= ['obsHistID', 'expMJD', \
	    'fieldRA','fieldDEC']
# identify column index of parameters to be retrieved
# + identify MJD and airmass columns
hecolist = []
for hh in headextra:
    hecolist = hecolist + [head.index(hh)]
    if hh == 'expMJD' :
        colmjd = head.index(hh)
    if hh == 'fieldRA' :
        colRA = head.index(hh)
    if hh == 'fieldDEC' :
        colDec = head.index(hh)
# 
#
#
# extract usefull data and format in single string
#      field data are stored in table :
#                 pplis
#     as single strings in the format :
#           parm_0 = val_0 $ parm_1 = val_1 $ ....
#      numerical values of mjd and and zenith angle
#      are stored in eponym tables
#

obshlines = fullhlines
#
print 'Succesfully read ', len(obshlines),' pointings in file : '\
      ,LSSTdatafile 
  
nobs = len(obshlines)
pplis = [[]]*nobs
mjdlis = [[]]*nobs
mjdlis[0] = MJDS
RAlis = [[]]*nobs
Declis = [[]]*nobs
zanglis = [[]]*nobs
for ipl in range(1,nobs):
    pl = obshlines[ipl][:-1]
    pli=''
    pline=pl.split('\t')
    for ic in hecolist:
        pli = pli + head[ic] + ' = ' + pline[ic] + ' $ '
        pplis[ipl] = pli
        mjdlis[ipl] = float(pline[colmjd])
        RAlis[ipl] = float(pline[colRA])
        Declis[ipl] = float(pline[colDec])
    nend = ipl+1
    #print nend
print 'pointing parms extracted for  ', nend,'  pointings'
# pplis contains one parameter string per LSST pointing
# mjdlis contains MJD times of LSST pointings 
#_______________________________________________________________________________
#=====================================================
# Search for atmospheric parameters at pointing times
#   derive zangle and azimuth from RA, Dec, MJD
#   derive aerosol parameters from list
#   format pointing sequence
#
for ipoint in range(1,nend) :
    Ra = RAlis[ipoint]*r2d
    Dec = Declis[ipoint]*r2d
    mjd = mjdlis[ipoint]
    azza = Astrotools.eq2loc([Ra,Dec,mjd])
    azim = azza[0] # field azimuth at mjd(degrees)
    zang = azza[1] # field zenith angle at mjd (degrees)
    idmjd = int((mjd-MJDS)/mjdstep)
    atmoparm = atpl[idmjd][0]
    vis0 = avpl[idmjd][0]
    vamp = avpl[idmjd][1] # origine des azimuth pour les aerosols
    vphi = avpl[idmjd][2]
    zaer11 = avpl[idmjd][3]
    zaer12 = avpl[idmjd][4]
    scale1 = avpl[idmjd][5]
    zaer21 = avpl[idmjd][6]
    zaer22 = avpl[idmjd][7]
    scale2 = avpl[idmjd][8]
    si2z = 1.
    if zang < 45. :
        si2z = sin(2 *zang * d2r)
    vis = vis0 + vamp * sin(azim * d2r - vphi)* si2z
    if vis < 4.0:
        vis = 4.0
    pplis[ipoint] = 'ID' + pplis[ipoint][12:]\
                + ' ZANGLE = % 8.3f' %zang\
                + ' $ '+ atmoparm \
                + ' VIS = % 8.1f'%vis\
                + ' $ ZAER11 = % 6.0f'%zaer11\
                + ' $ ZAER12 = % 6.0f'%zaer12\
                + ' $ SCALE1 = % 6.0f'%scale1\
                + ' $ ZAER21 = % 6.0f'%zaer21\
                + ' $ ZAER22 = % 6.0f'%zaer22\
                + ' $ SCALE2 = % 6.0f'%scale2 + ' $ \n'

######################################################
#Write pointing sequence with atmosphere parameters
#  
sequence_Parmfile = mlsstpath +'Runparm/'\
                    + expname + '_parmlist.dat'
seqparmfil = open(sequence_Parmfile,'w')
for ic in range(1,nend):
    seqparmfil.write(pplis[ic])
seqparmfil.close()
print nend,'parameter lines  written in batch ',expname
print '       from ',pplis[1][:20], '  to ',pplis[ipoint][:20]
