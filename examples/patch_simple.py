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
import lsst.sims.selfcal.analysis.catalog_plots_utils as cpu
import lsst.sims.selfcal.analysis.catalog_plots as cp
import os.path
import subprocess as sb

pyl.rcParams['font.size']=16


# true patch data file (usually patchdata.dat)
# NOTE THAT 'patchdata.dat.s' is the SORTED (by patch ID) version of the patchdata.dat file

truepatchfile = "patchdata.dat.s"
if not os.path.isfile(truepatchfile):
    sor=sb.Popen('sort -nk 2 patchdata.dat > patchdata.dat.s', shell=True).wait()


# best fit patch data file (usually test_bestfit_Patch.dat)
if len(sys.argv)<2:
    bestfitpatchfile = "test_bestfit_Patch.dat"
else:
    bestfitpatchfile = sys.argv[1]

iteration_number = -99
delta_chisq = -99 

if len(sys.argv)>3:
    iteration_number = int(sys.argv[2])
    delta_chisq = float(sys.argv[3])

print "Using %s as the true patch data file, and %s as the bestfit patch data file" \
      %(truepatchfile, bestfitpatchfile)

print "also iteration number %d and delta_chisq %g" %(iteration_number, delta_chisq)

# read the data files from disk
patchdat = cp.read_true_patchfile(truepatchfile)
patchcal = cp.read_bestfit_patchfile(bestfitpatchfile)

# match the true and best-fit patch
patchdat = cp.match_patch_calvstrue(patchdat, patchcal, sorted=True)

# adjust the calibrated patch so they have approximately the right overall zeropoint
patchdat = cp.adjust_patch_calvstrue(patchdat)
if n.size(n.where(patchdat['patchid'] == 0)) !=0:
    good = n.where(patchdat['patchid'] == 0)
    print 'Zeropoint from flux calibration stars = %f'%patchdat['magdiff'][good]
# plot the difference between the bestfit and true patch magnitudes
maglim = None   # or set magnitude limit ranges 
# or set magnitude limit ranges

#use a 90% cutoff to clip outliers that skew the stretch
ord=n.abs(patchdat['magdiff']).argsort()
ord=ord[0:len(ord)*.9]

rms = patchdat['magdiff'][ord].std()

print 'zp - fit RMS = %f, zp - fit (clipped) = %f'%(patchdat['magdiff'].std(), rms)

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

cp.plot_patch_dmagcaltrue(patchdat, maglim=maglim, etitle=etitle, z_method=n.mean)
savefig=True
figformat = "png"
if iteration_number > 0 :
    figname = "P%ddmag." %(iteration_number)
else:
    figname = "Pdmag."
if savefig:
    pyl.savefig(figname + figformat, format=figformat)

maglim2 = [0, round(maglim[1]/3*1000)/1000.0]
#cp.plot_patch_dmagcaltrue(patchdat, maglim=maglim2, etitle=" RMS "+etitle, z_method=n.std)
#if iteration_number > 0 :
#    figname = "P%ddmag_rms." %(iteration_number)
#else:
#    figname = "Pdmag_rms."
#if savefig:
#    pyl.savefig(figname + figformat, format=figformat)

def delta(z):
    d = n.max(z) - n.min(z)
    return d
maglim2 = [0, maglim[1]/2]
#cp.plot_patch_dmagcaltrue(patchdat, maglim=maglim2, etitle=" DELTA "+etitle, z_method=delta)
#if iteration_number > 0 :
#    figname = "P%ddmag_delta." %(iteration_number)
#else:
#    figname = "Pdmag_delta."
#if savefig:
#    pyl.savefig(figname + figformat, format=figformat)


# plot the original offsets applied too
#maglim_zpin = [0, 1.5]
#cp.plot_patch_zpoffset(patchdat, maglim=maglim_zpin, etitle=etitle, z_method=n.mean)
#if iteration_number > 0 :
#    figname = "P%dzp." %(iteration_number)
#else:
#    figname = "Pzp."
#if savefig:
#    pyl.savefig(figname + figformat, format=figformat)

#cp.plot_patch_zpoffset(patchdat, maglim=maglim_zpin, etitle=" MAX " + etitle, z_method=n.max)
#if iteration_number > 0 :
#    figname = "P%dzp_min." %(iteration_number)
#else:
#    figname = "Pzp_min."
#if savefig:
#    pyl.savefig(figname + figformat, format=figformat)

#cp.plot_patch_zpoffset(patchdat, maglim=maglim2, etitle=" RMS " + etitle, z_method=n.std)
#if iteration_number > 0 :
#    figname = "P%dzp_rms." %(iteration_number)
#else:
#    figname = "Pzp_rms."
#if savefig:
#    pyl.savefig(figname + figformat, format=figformat)


outlier_clip=n.where(n.abs(patchdat['magdiff']) < .05)
rms = patchdat['magdiff'][outlier_clip].std()


# also plot the histogram of the difference in true and bestfit magnitudes
if maglim==None:
    histrange= None
else:
    histrange = [maglim[0]*2, maglim[1]*2]  # or set range for histogram
nbins = 100       # or set number of bins for histogram
cp.histogram_patch_dmagcaltrue(patchdat, nbins=nbins, range=histrange, etitle=etitle+"RMS %4.2f mmag"%(1000.*rms))
if iteration_number > 0 :
    figname = "P%ddmag_hist." %(iteration_number)
else:
    figname = "Pdmag_hist."
if savefig:
    pyl.savefig(figname + figformat, format=figformat)

#if maglim_zpin==None:
#    histrange_zpin=None
#else:
#    histrange_zpin = [maglim_zpin[0], maglim_zpin[1]]
#cp.histogram_patch_zpoffset(patchdat, nbins=nbins, range=histrange_zpin, etitle=etitle)
#if iteration_number > 0 :
#    figname = "P%dzp_hist." %(iteration_number)
#else:
#    figname = "Pzp_hist."
#if savefig:
#    pyl.savefig(figname + figformat, format=figformat)
    
#pyl.show()

#plot up the number of stars per patch
pyl.figure()
bins=n.arange(n.median(patchdat['nstars'])*3)
num, b,p = pyl.hist(patchdat['nstars'], bins=bins, normed=1)
pyl.title('Number of stars per patch, median=%.0f'%(n.median(patchdat['nstars'])))
pyl.xlabel('Number of stars')
#pyl.ylabel('')
pyl.savefig('Pnstars_hist.png', format='png')


#cp.plot_patch_nstars(patchdat, maglim=[0,round(n.median(patchdat['nstars']))*3])
#pyl.savefig('Pnstars.png',format='png')

#calculate an RMS for each sub-patch
#n_sub=n.max(patchdat['subpatch'])+1
#nside=n.ceil(n_sub**0.5)
#n_sub=nside**2
#subs_rms=n.zeros(n_sub+0.)
#for i in n.arange(n_sub):
#    magdiff=patchdat['magdiff'][n.where(patchdat['subpatch'] == i)]
    #do a 90th percentile clip to get rid of outliers
#    if n.size(magdiff) > 5:
#        ord=n.abs(magdiff).argsort()
#        ord=ord[0:len(ord)*.9]
#        subs_rms[i] = magdiff[ord].std()
#now need to assign radius to each subpatch, and then color code it.
#nside=n_sub**0.5
#image=n.zeros((nside,nside))
#for i in n.arange(nside):
#    for j in n.arange(nside):
#        image[i,j]=subs_rms[i*nside+j]*1000

#pyl.figure()
#pyl.imshow(image, extent=(0,nside,0,nside), interpolation='nearest')
#pyl.title('Subpatch RMS')
#cb = pyl.colorbar()
#cb.set_label('mmag')
#pyl.savefig('Psubpatch.png',type='png')
