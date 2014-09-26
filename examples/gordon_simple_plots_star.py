#!/usr/bin/env python

#####
#  Lynne Jones, ljones@astro.washington.edu.
#  svn version info : $Id: simple_plots.py 16346 2010-07-15 19:59:35Z ljones $
#
# This python script will read the bestfit and true star data files, and create plots
# of their differences. 
#
#####

import sys

import numpy as n
import matplotlib
matplotlib.use('Agg')
import pylab as pyl

import lsst.catalog_plots_utils as cpu
import lsst.catalog_plots as cp

pyl.rcParams['font.size']=16
#pyl.rcParams['axes.labelsize']='large'
#pyl.rcParams['xtick.labelsize']='large'
#pyl.rcParams['ytick.labelsize']='large'


deg2rad = n.pi/180.0



# true star data file (usually stardata.dat)
truestarfile = "stardata.dat"
# best fit star data file (usually test_bestfit_Star.dat)
if len(sys.argv)<2:
    bestfitstarfile = "test_bestfit_Star.dat"
else:
    bestfitstarfile = sys.argv[1]

iteration_number = -99
delta_chisq = -99 

if len(sys.argv)>3:
    iteration_number = int(sys.argv[2])
    delta_chisq = float(sys.argv[3])

print "Using %s as the true star data file, and %s as the bestfit star data file" \
      %(truestarfile, bestfitstarfile)

print "also iteration number %d and delta_chisq %g" %(iteration_number, delta_chisq)

# read the data files from disk
stardat = cp.read_true_stars(truestarfile) #keys are ('id', 'dbid', 'ra', 'dec', 'magtrue', 'color')
starcal = cp.read_bestfit_stars(bestfitstarfile) #kets are ('id', 'magcal')

# match the true and best-fit stars
stardat = cp.match_stars_calvstrue(stardat, starcal)

# adjust the calibrated stars so they have approximately the right overall zeropoint
stardat = cp.adjust_stars_calvstrue(stardat)

# plot the difference between the bestfit and true star magnitudes
maglim = None
# or set magnitude limit ranges
outlier_clip=n.where(n.abs(stardat['magdiff']) < .05)
rms = stardat['magdiff'][outlier_clip].std()
print 'true-bestfit RMS = %f, true-bestfit (clipped) = %f'%(stardat['magdiff'].std(), rms)
lim = (rms*2.5)
if lim < 0.001:
    lim = round(lim * 10000) / 10000.0
elif lim < 0.01:
    lim = round(lim * 1000) / 1000.0
print "Using magnitude limits for plots of %f, %f" %(-lim, lim)
maglim= [-lim, lim]

# set extra title words for extra description
if iteration_number > 0:
    etitle = " Iter: %d  Delta_Chi_sq: %.2g" %(iteration_number, delta_chisq) 
else:
    etitle = ""
   
cp.plot_stars_dmagcaltrue(stardat, maglim=maglim,
                          etitle=etitle+"RMS %f mmag"%(1000.*stardat['magdiff'].std()), z_method=n.mean )

savefig=True
figformat = "png"
if iteration_number > 0 :
    figname = "S%ddmag." %(iteration_number)
else:
    figname = "Sdmag."
if savefig:
    pyl.savefig(figname + figformat, format=figformat)

if maglim==None:
    maglim = None
else:
    maglim2 = [0, maglim[1]/5.] #was [0, maglim[1]*1.5]
cp.plot_stars_dmagcaltrue(stardat, maglim=maglim2, etitle=" RMS "+etitle, z_method=n.std)

if iteration_number > 0 :
    figname = "S%ddmag_rms." %(iteration_number)
else:
    figname = "Sdmag_rms."
if savefig:
    pyl.savefig(figname + figformat, format=figformat)

# also plot the histogram of the difference in true and bestfit magnitudes
if maglim == None:
    histrange = None
else:
    histrange = [maglim[0]*2, maglim[1]*2]  # or set range for histogram
nbins = 100       # or set number of bins for histogram
cp.histogram_stars_dmagcaltrue(stardat, nbins=nbins, range=histrange, etitle=etitle+"RMS %f mmag"%(1000.*stardat['magdiff'].std()))

if iteration_number > 0 :
    figname = "S%ddmag_hist." %(iteration_number)
else:
    figname = "Sdmag_hist."
if savefig:
    pyl.savefig(figname + figformat, format=figformat)


cpu.plot_density(stardat['x'],stardat['y'],z=stardat['color'],radecgrid=True)
pyl.title('Average $g-i$')
if savefig:
    pyl.savefig('Starcolor.png', format='png')


cpu.plot_density(stardat['x'],stardat['y'],z=stardat['color'],z_method=n.max, radecgrid=True)
pyl.title('max $g-i$')
if savefig:
    pyl.savefig('Starcolormax.png', format='png')

cpu.plot_density(stardat['x'],stardat['y'],z=stardat['color'],z_method=n.min,radecgrid=True)
pyl.title('min $g-i$')
if savefig:
    pyl.savefig('Starcolormin.png', format='png')

cpu.plot_density(stardat['x'],stardat['y'],z=stardat['color'],z_method=n.std,radecgrid=True)
pyl.title('$\sigma(g-i)$')
if savefig:
    pyl.savefig('Starcolorstd.png', format='png')
pyl.figure()
ack1,ack2,ack3 = pyl.hist(stardat['color'], bins=100)
pyl.xlabel('$(g-i)$')
pyl.ylabel('Number of stars')
pyl.title('Color Distribution of Calibration Stars')
pyl.savefig('Starcolorhist.png', format='png')

