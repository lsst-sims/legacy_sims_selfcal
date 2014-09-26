#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#"""                                                                   """
#"""                        Dmag_sequenceV5-multiref.py :                         """
#"""                                                                   """
#          Compute magnitude variations due to atmosphere
#   this is same as V4 but the grid of the reference model
#    at various airmasses is equi_spaced in arimass (step0.01)
# instead of in zangle (degree)
#--------------------------------------------------------------------------
# Read :
#       - filter responses
#       - reference transmission spectrum given at each degree of zangle
#       - pointing identification and airmass of each pointing
#       - Modtran atmosphere spectra for the sequence of pointings
#   
# Compute : magnitudes differences (pointing-reference (at same airmass))
#
# Write in _dmag.dat file
#=========================================================================
import numpy
from numpy import *
import Astrotools


localpath =  '/groups/LSST/creze/Calsimwork/'
mlsstpath = localpath +'mkatmos/SpData/'

expname = 'Opsim3.61.02sel'

refname = 'Ref400'
specpath = mlsstpath + expname + '/'
refpath  = mlsstpath + refname + '/'
runpath = mlsstpath + 'Runparm/'

phopath = localpath + 'mkatmos/Parmfiles/'
phofile = 'LSSTfiltersonMDTRbins.prn'



nsbin = 5212
#=========================================================================
#
#           FILTER DATA ACQUISITION (return as array LSfi)
#
# ------------------------------------------------------------------------
#
finames = '      filter           :       u       |        g       |        r       |'\
          + '       i       |        z       |        y \n'
fipho = phopath + phofile
LSSTfilf = open(fipho,'r')
LSfi=zeros((nsbin,6))
#
il=-1
for lf in LSSTfilf:
    il=il+1
    lulu = lf.split()
    lili = [float(lulu[1]), float(lulu[2]),\
            float(lulu[3]),float(lulu[4]),\
            float(lulu[5]), float(lulu[6])]
    LSfi[il,:]= array(lili)
                                                                       
LSSTfilf.close()
#
#                 END FILTER ACQUISITION
#
#======================================================================
#
#              TRANSMISSION SPECTRUM DATA ACQUISITION
#
#----------------------------------------------------------------------
#
#----------------
# Read experience sequence description ( return as "comment" )
#----------------
#
Descriptfile = runpath + expname+'_descript.txt'
descript = open(Descriptfile,'r')
tdesc = []
for desc in descript:
    tdesc = tdesc + [desc]
descript.close() # full description file is now in list tdesc
idc1 = tdesc.index('SEQUENCE NAME\n')+ 1
idc2 = tdesc.index('SEQUENCE DESCRIPTION\n')+1
comment = [tdesc[idc1],tdesc[idc2],tdesc[-2],tdesc[-1]] #  text desription is now in list comment
#
#---------
# Read Sequence Pointing identifiers and airmass
#            from runparm ..._parmlist.dat)
#---------
#
Seqfile= runpath + expname + '_parmlist.dat'
sfi = open(Seqfile,'r')
 # initiate observation identification list
obsidl = ['  Reference ']
airmass = [0.] # initiate airmass list
airmass = airmass + [1.]
zidx = []
zidx = zidx + [0]
ipp = 0
# search 
for parml in sfi:
	ipp= ipp+1
	if ipp == 1 :
		parl = parml.split('$')
		ident = parl[0]
		obsidl = obsidl+[ident]
		for ii in range (1,len(parl)) :
			pa = parl[ii].split('=')
			if pa[0].strip() == 'ZANGLE':
				iam = ii
				zangle = numpy.float(pa[1])
				am = Astrotools.mdtram(zangle)
				airmass = airmass+[am]
				zidx = zidx + [numpy.int((am-1)*100.)]
	else :
		parl = parml.split('$')
		obsidl =obsidl + [ parl[0]]
		pa = parl[iam].split('=')
		zangle = numpy.float(pa[1])
		am = Astrotools.mdtram(zangle)
		airmass = airmass+[am]
		zidx = zidx + [numpy.int((am-1)*100.)]
#---------
# Read spectrum data for Reference
#---------
#
firef = refpath + refname + 'csp1.dat'
refsp = open(firef,'r')
refspli = list(refsp)
refsp.close()
alj = []
for lila in refspli:
    lila = lila.split('$')
    for ix in range(1,len(lila)):# lila[0] contains wavelengthes
        lila[ix]=float(lila[ix])
    alj = alj +[lila[1:]]
Reference = array(alj)
#---------
# Read spectrum data for sequence (from c(ompact)sp(ectra)  files)
#---------
#

ficont = specpath + expname +'csp_control.dat'
fic = open(ficont,'r')
spfides = list(fic)
for desc in spfides :
    descl = desc.split('=')
    if descl[0].strip() == 'block number':
        blknum=int(descl[1])
#========================================================
#  initialize result file
#
resufil= open(specpath+expname+'_dmag.dat','w')
for com in comment:
    resufil.write(com)   # write comment from descript file
resufil.write(finames)   # write filter names 
#========================================================
#  process data by block files
#
for ibl in range(0,blknum):
    ns00 = ibl*1000
    bloname = specpath + expname + 'csp' + str(ibl+1)+'.dat'
    trsp = open(bloname,'r')
    lignes = list(trsp)
    if len(lignes) <> nsbin :
        print '***** number of bins mismatch*****'
    # convert charline data spectrum to numerical array spectrum
    ali = []
    for lila in lignes:
        lila = lila.split('$')
        for ix in range(0,len(lila)):# lila[0] contains wavelengthes
            lila[ix]=float(lila[ix])
        ali = ali +[lila[1:]]
    Lamspec = array(ali)
    #
    #               END TRANSMISSION SPECTRA ACQUISITION
    #
    #======================================================================
    #
    #                 COMPUTE MAGNITUDE DIFFERENCES 
    #
    #----------------------------------------------------------------------
    namelist = []
    #   Compute magnitude deviations with respect to standard atmos at
    #    same zangle
    #  nspec = 10
    nspec = shape(Lamspec)[1]
    for isp in range(0,nspec):
        jsp = ns00 + isp
        if isp==0 :
            jsp = isp
        tf  = zeros((nsbin,6))# init. transm. through filters + current Atm.
        t0f = zeros((nsbin,6))# init. transm. through filters + standard Atm.
        Ref = Reference[:,zidx[jsp]] # select reference for current airmass
        ef= zeros(6) # initialize total energy through fi +current
        e0f = zeros(6)# initialize total energy through fi+standard
        for i in range(0,6):
            tf[:,i] = Lamspec[:,isp]*LSfi[:,i]# energy in bins through current A.
            t0f[:,i] = Ref * LSfi[:,i]# energy in bins through standard A.
        ef = tf.sum(axis=0)# energy in filter through atmosphere
        e0f = t0f.sum(axis=0)# energy in filter through atmosphere
        dmag = -2.5 * log10 (ef/e0f) #magnitude variation current/standard
        obid = obsidl[isp+ns00]
        if isp == 0 :
            obid = ' Reference'
        resline1 = ' '+obid+' :\n'
#
        for i in range(0,6):
            resline1 = resline1[:-1] +'%8.4f'%dmag[i]+'    \n'
        resufil.write(resline1)
resufil.close()



