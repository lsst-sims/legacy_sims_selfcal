#!/usr/bin/env python

####
# python script to make various metric plots based on the
#output of self-calibration run on HEALpixels
###

import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as pyl
from scipy import stats
import lsst.sims.selfcal.analysis.catalog_plots_utils as cpu
import lsst.sims.selfcal.analysis.catalog_plots as cp
import lsst.sims.selfcal.generation.split_obs as so
import lsst.sims.selfcal.analysis.useful_input as ui
import healpy as hp

pyl.rcParams['font.size']=16
deg2rad = np.pi/180.0


def dictsort(indict,key):
    """for a dictionary of numpy arrays.  Normally, one would just use record arrays and np.sort(array,order=key) """
    order = np.argsort(indict[key])
    for key in indict.keys():
        indict[key] = indict[key][order]
    return indict



#read in patch data, match it up, and make some plots.
def patchPlots(truepatchfile='patchdata.dat', bestfitpatchfile = "test_bestfit_Patch.dat"):
    patchdat = dictsort(cp.read_true_patchfile(truepatchfile), 'patchid')
    patchcal = cp.read_bestfit_patchfile(bestfitpatchfile)
    patchdat = cp.match_patch_calvstrue(patchdat, patchcal, sorted=True)
    #Remove the floating zeropoint
    patchdat = cp.adjust_patch_calvstrue(patchdat)
    output = {}
    output['N Patches Fit'] = np.size(patchdat['patchid'])
    output['N visits'] = np.size(np.unique(patchdat['visitid']))
    #use a 90% cutoff to clip outliers that skew the stretch
    ord=np.abs(patchdat['magdiff']).argsort()
    ord=ord[0:len(ord)*.9]

    rms = patchdat['magdiff'][ord].std()
    output['Residuals RMS (mmag)'] = rms*1e3
    output['IQR_sigma (mmag)'] = so.robustSigma( patchdat['magdiff'])*1e3
    print 'zp - fit RMS = %f, zp - fit (clipped) = %f'%(patchdat['magdiff'].std(), rms)
    lim = (rms*2.5)
    if lim < 0.001:
        lim = round(lim * 10000) / 10000.0
    elif lim < 0.01:
        lim = round(lim * 1000) / 1000.0
    print "Using magnitude limits for plots of %f, %f" %(-lim, lim)
    maglim= [-lim, lim]
    #residual map
    etitle=""
    figname = "Pdmag."
    figformat='png'
    cp.plot_patch_dmagcaltrue(patchdat, maglim=maglim, etitle=etitle, z_method=np.mean)
    pyl.savefig(figname + figformat, format=figformat)


    if maglim==None:
        histrange= None
    else:
        histrange = [maglim[0]*2, maglim[1]*2]  # or set range for histogram
    nbins = 100       # or set number of bins for histogram
    cp.histogram_patch_dmagcaltrue(patchdat, nbins=nbins, range=histrange, etitle=etitle+"RMS %4.2f mmag"%(1000.*rms))
    figname = "Pdmag_hist."
    pyl.savefig(figname + figformat, format=figformat)


    pyl.figure()
    bins=np.arange(np.median(patchdat['nstars'])*3)
    num, b,p = pyl.hist(patchdat['nstars'], bins=bins)#, normed=1)
    pyl.title('Number of stars per patch, median=%.0f'%(np.median(patchdat['nstars'])))
    pyl.xlabel('Number of stars')
    pyl.savefig('Pnstars_hist.png', format='png')
    output['ave starsPerPatch'] = np.median(patchdat['nstars'])


    pyl.figure()
    num, b,p = pyl.hist(patchdat['dmagcal']-np.min(patchdat['dmagcal']), bins=100)#, normed=1)
    pyl.xlabel('Patch Zeropoint')
    pyl.ylabel('Number of Patches')
    output['Photometric Fraction'] = float(np.size(np.where(patchdat['dmagcal']-np.min(patchdat['dmagcal']) < .01)))/np.size(patchdat['dmagcal'])
    pyl.title('Extinction Distribution')
    pyl.text(.4,.8, 'Photometric Fraction = %3.3f'%output['Photometric Fraction'], transform=pyl.gca().transAxes)
    pyl.savefig('Pzphist.png', format='png')
    


    return output


def starPlots(truestarfile = "stardata.dat", bestfitstarfile = "test_bestfit_Star.dat"):

    #assume both the files are sorted by star id.
    stardat = cp.read_true_stars(truestarfile) #keys are ('id', 'dbid', 'ra', 'dec', 'magtrue', 'color')
    starcal = ui.readDatafile(bestfitstarfile, keys=('id','magcal','rms','iqr_sigma','nobs'))
    HP = True
    if np.size(starcal['id']) == 0:
        starcal = ui.readDatafile(bestfitstarfile, keys=('id','magcal'))
        HP = False                     
    fitted_stars = np.in1d(stardat['id'], starcal['id'])
    for key in stardat.keys():
        stardat[key] = stardat[key][fitted_stars]
    for key in starcal.keys():
        stardat[key] = starcal[key]
    stardat['magdiff'] = stardat['magcal'] - stardat['magtrue']
    stardat = cp.adjust_stars_calvstrue(stardat)

    starResults = {}
    HPResults = {}
    starResults['N stars Fit'] = np.size(stardat['magdiff'])

    outlier_clip=np.where(np.abs(stardat['magdiff']) < .05)
    rms = stardat['magdiff'][outlier_clip].std()
    print 'true-bestfit RMS = %f, true-bestfit (clipped) = %f'%(stardat['magdiff'].std(), rms)
    lim = rms*2.5
    if lim < 0.001:
        lim = round(lim * 10000) / 10000.0
    elif lim < 0.01:
        lim = round(lim * 1000) / 1000.0
    print "Using magnitude limits for plots of %f, %f" %(-lim, lim)
    maglim= [-lim, lim]
    
    starResults['Residuals RMS (mmag)'] = stardat['magdiff'].std()*1e3
    starResults['IQR_sigma (mmag)'] = so.robustSigma(stardat['magdiff'])*1e3
    
    #star residual map
    stardat['x'], stardat['y'] = cpu.hammer_project_toxy(stardat['ra']*deg2rad,
                                                         stardat['dec']*deg2rad)
    gridsize = cpu.calc_gridsize(stardat, binsize='patch')
    etitle = " RMS = %4.1f mmag"%(rms*1000)
    cpu.plot_density(stardat['x'], stardat['y'], stardat['magdiff']*1000, z_method=np.mean,
                     zlimits=np.array(maglim)*1000, gridsize=gridsize, radecgrid=True, newfig=True, cb_label='mmag')
    figtitle = "Stars dMag(true-bestfit)"
    pyl.title(figtitle, fontsize='medium')
    pyl.savefig('Sdmag.png',type='png')

    
    maglim2 = [0, maglim[1]/2.]
    cp.plot_stars_dmagcaltrue(stardat, maglim=maglim2, etitle=" RMS "+etitle, z_method=np.std)
    pyl.savefig('Sdmag_rms.png', format='png')

    histrange = [maglim[0]*2.5, maglim[1]*2.5]
    nbins = 100
    cp.histogram_stars_dmagcaltrue(stardat, nbins=nbins, range=histrange, etitle=etitle)
    pyl.savefig('Sdamg_hist.png', format='png')


   
  
    #color distribution of stars
    pyl.figure()
    pyl.hist(stardat['color'], bins=100)
    pyl.xlabel(r'$g-i$')
    pyl.ylabel('# of stars')
    pyl.title('Color Distribution')
    pyl.savefig('Scolordist.png',type='png')


    #accuracy v color
    pyl.figure()
    #pyl.scatter(stardat['color'], stardat['magdiff']*1000., c=stardat['magtrue'], alpha=0.1,edgecolors=None)
    pyl.hexbin(stardat['color'], stardat['magdiff']*1000., bins='log')
    cb=pyl.colorbar()
    cb.ax.set_ylabel(r'$\log{\rm{N}}$')
    pyl.xlabel(r'$g-i$')
    pyl.ylabel('Bestfit-True (mmag)')
    pyl.savefig('Saccuracyvcolor.png',type='png')

    
    #accuracy v mag
    pyl.figure()
    #pyl.scatter(stardat['magtrue'], stardat['magdiff']*1000., c=stardat['color'], alpha=0.1,edgecolors=None)
    pyl.hexbin(stardat['magtrue'], stardat['magdiff']*1000., bins='log')
    cb=pyl.colorbar()
    cb.ax.set_ylabel(r'$\log{\rm{N}}$')
    pyl.xlabel('True Mag')
    pyl.ylabel('Bestfit-True (mmag)')
    pyl.savefig('Saccuracyvmag.png',type='png')

    if HP:

        zlim = [np.min(stardat['nobs']), (np.median(stardat['nobs'])-np.min(stardat['nobs']))+np.median(stardat['nobs'])]
        cpu.plot_density(stardat['x'], stardat['y'],stardat['nobs'], z_method=np.mean,
                         zlimits=zlim, gridsize=gridsize, radecgrid=True, newfig=True, cb_label='N observations')
        pyl.title('Visit Density')
        pyl.savefig('Snobs.png', format='png')

        pyl.figure()
        #I should just figure out from the data what nsides is!
        hpid = hp.ang2pix(16, (stardat['dec']+90.)*deg2rad, stardat['ra']*deg2rad)
        cpu.plot_density(stardat['x'], stardat['y'],np.mod(hpid,250), z_method=np.median,
                          gridsize=gridsize, radecgrid=True, newfig=True, cb_label='HEALpix ID mod 250')
        pyl.title("HEALpix Map")
        pyl.savefig('HPid.png')
        HPResults['number of HEALpix'] = np.size(np.unique(hpid))
        
        pyl.figure()
        num, b, p = pyl.hist(stardat['nobs'], bins=100, range=zlim, histtype='bar')
        pyl.xlabel("Number of Observations per Star")
        pyl.title('Median = %d'%np.median(stardat['nobs']))
        pyl.savefig('Snobs_hist.png',type='png')
        starResults['median repeat obs'] = np.median(stardat['nobs'])

        pyl.figure()
        x = np.sort(stardat['iqr_sigma'])*1e3
        y = np.arange(np.size(x), dtype='float')
        y = y/np.max(y)
        per50 = np.max(np.where(y <= .5))
        per90 = np.max(np.where(y <= .9))
        pyl.plot(x,y, 'k--', label='IQR, 50th,90th %2.1f, %2.1f'%(x[per50],x[per90]))
        pyl.plot([0,x[per50]],[y[per50],y[per50]],'b--')
        pyl.plot([x[per50],x[per50]],[0,y[per50]],'b--')
        pyl.plot([0,x[per90]],[y[per90],y[per90]],'b--')
        pyl.plot([x[per90],x[per90]],[0,y[per90]],'b--')
        #pyl.title('Cumulative Distribution of Stellar Repeatability')
        #pyl.xlabel('Repeatability IQR RMS (mmag)')
        #pyl.ylabel('Cumulative Fraction')
        #pyl.savefig('Srepeat_IQR_cumulative.png',type='png')
        pyl.legend()
        pyl.savefig('Srepeat_cumulative.png',type='png')

       #repeatability from IQR
        rs = so.robustSigma(stardat['iqr_sigma'])
        med = np.median(stardat['iqr_sigma'])
        maglim = np.around(np.array([med-3*rs, med+3.*rs])*1000, decimals=2)/1000
        cpu.plot_density(stardat['x'], stardat['y'], stardat['iqr_sigma']*1000, z_method=np.mean,
                         zlimits=maglim*1000, gridsize=gridsize, radecgrid=True, newfig=True, cb_label='mmag')
        pyl.title('Repeatability (IQR)')
        pyl.savefig('Srepeat_IQR.png',type='png')


        pyl.figure()
        num, b, p = pyl.hist(stardat['iqr_sigma']*1000, bins=100, range=maglim*1000, histtype='bar')
        pyl.xlabel('RMS per star (mmag)')
        pyl.title('Repeatability IQR, Median = %4.2f mmag'%np.median(stardat['iqr_sigma']*1000))
        pyl.savefig('Srepeat_IQR_hist.png',type='png')
        starResults['Median Repeatability (IQR) (mmag)'] = np.median(med)*1e3


        pyl.figure()
        good = np.where(stardat['magtrue'] < 18)
        num, b, p = pyl.hist(stardat['iqr_sigma'][good]*1000, bins=100, range=maglim*1000, histtype='bar')
        pyl.xlabel('RMS per star, m < 18 (mmag)')
        pyl.title('Repeatability IQR, Median = %4.2f mmag'%np.median(stardat['iqr_sigma'][good]*1000))
        pyl.savefig('Srepeat_IQR_bright_hist.png',type='png')
        starResults['Median Repeatability Bright (IQR) (mmag)'] = np.median(stardat['iqr_sigma'][good])*1e3
     

        #repeatability v color
        pyl.figure()
        #pyl.plot(stardat['color'], stardat['rms']*1000., 'ko', alpha=.1)
        pyl.hexbin(stardat['color'], stardat['iqr_sigma']*1000., bins='log')
        cb = pyl.colorbar()
        cb.set_label(r'$\log{\rm{N}}$')
        pyl.ylim([0,40])
        pyl.xlabel(r'$g-i$')
        pyl.ylabel(r'Repeatability $\sigma_{\rm{IQR}}$ (mmag)')
        pyl.savefig('Srepeatvcolor.png',type='png')

        #repeatability v mag
        pyl.figure()
        #pyl.plot(stardat['magtrue'], stardat['rms']*1000., 'ko', alpha=.1)
        pyl.hexbin(stardat['magtrue'], stardat['iqr_sigma']*1000.,bins='log')
        cb = pyl.colorbar()
        cb.set_label(r'$\log{\rm{N}}$')
        pyl.ylim([0,40])
        pyl.xlabel('True Mag')
        pyl.ylabel(r'Repeatability $\sigma_{\rm{IQR}}$ (mmag)')
        pyl.savefig('Srepeatvmag.png',type='png')


        #repeatability
        rs = so.robustSigma(stardat['rms'])
        med = np.median(stardat['rms'])
        maglim = np.around(np.array([med-3*rs, med+3.*rs])*1000, decimals=2)/1000
        cpu.plot_density(stardat['x'], stardat['y'], stardat['rms']*1000, z_method=np.mean,
                         zlimits=maglim*1000, gridsize=gridsize, radecgrid=True, newfig=True, cb_label='mmag')
        pyl.title('Repeatability')
        pyl.savefig('Srepeat.png',type='png')

        pyl.figure()
        num, b, p = pyl.hist(stardat['rms']*1000, bins=100, range=maglim*1000, histtype='bar')
        pyl.xlabel('RMS per star (mmag)')
        pyl.title('Repeatability, Median = %4.2f mmag'%np.median(stardat['rms']*1000))
        pyl.savefig('Srepeat_hist.png',type='png')
        starResults['Median Repeatability (mmag)'] = np.median(med)*1e3

        pyl.figure()
        x = np.sort(stardat['rms'])*1e3
        y = np.arange(np.size(x), dtype='float')
        y = y/np.max(y)
        per50 = np.max(np.where(y <= .5))
        per90 = np.max(np.where(y <= .9))
        pyl.plot(x,y, 'k', label='RMS, %2.1f, %2.1f'%(x[per50],x[per90]))
        pyl.plot([0,x[per50]],[y[per50],y[per50]],'k')
        pyl.plot([x[per50],x[per50]],[0,y[per50]],'k')
        pyl.plot([0,x[per90]],[y[per90],y[per90]],'k')
        pyl.plot([x[per90],x[per90]],[0,y[per90]],'k')
        pyl.title('Cumulative Distribution of Stellar Repeatability')
        pyl.xlabel('Repeatability RMS (mmag)')
        pyl.ylabel('Cumulative Fraction')
        pyl.savefig('Srepeat_cumulative.png',type='png')


        #need to make a spatial uniformity plot since we can now have varying densities of stars.  
        nside = 64
        stardat['hpid'] = hp.ang2pix(nside, (stardat['dec']+90.)*deg2rad, stardat['ra']*deg2rad)
        uhpid = np.unique(stardat['hpid'])
        hp_mean = uhpid*0.
        hp_std = uhpid*0.
        hp_dec, hp_ra = hp.pix2ang(nside,uhpid)
        #hp_ra = hp_ra
        hp_dec = hp_dec-90.*deg2rad
        hp_x, hp_y = cpu.hammer_project_toxy(hp_ra, hp_dec)
        stardat = dictsort(stardat,'hpid')

        left = np.searchsorted(stardat['hpid'], uhpid)
        right = np.searchsorted(stardat['hpid'], uhpid,side='right')

        for i in np.arange(left.size):
            hp_mean[i] = np.mean(stardat['magdiff'][left[i]:right[i]])
            hp_std[i] = np.std(stardat['magdiff'][left[i]:right[i]])

        cpu.plot_density(hp_x, hp_y, hp_mean*1000, z_method=np.mean,
                         zlimits=None, gridsize=gridsize, radecgrid=True, newfig=True, cb_label='mmag')
        HPResults['Spatial Uniformity RMS (mmag)'] = np.std(hp_mean)*1e3
        HPResults['Spatial Uniformity IQR RMS (mmag)'] = so.robustSigma(hp_mean)*1e3

        pyl.title('mean(True-Bestfit), binned by HEALpix')
        pyl.savefig('HPresid.png',type='png')

        cpu.plot_density(hp_x, hp_y, hp_std*1000, z_method=np.mean,
                         zlimits=None, gridsize=gridsize, radecgrid=True, newfig=True, cb_label='mmag')

        pyl.title('std(True-Bestfit), binned by HEALpix')
        pyl.savefig('HPstd.png',type='png')

        pyl.figure()
        num, b, p = pyl.hist( hp_mean*1000, bins=100, range=None, histtype='bar')
        pyl.xlabel("HEALpix(True-Bestfit) (mmag)")
        pyl.title('Spatial Uniformity, RMS=%4.2f'%HPResults['Spatial Uniformity RMS (mmag)'])
        pyl.savefig('HPmean_hist.png',type='png')

        pyl.figure()
        num, b, p = pyl.hist( hp_std*1000, bins=100, range=None, histtype='bar')
        pyl.xlabel("std(True-Bestfit) per HEALpix (mmag)")
        pyl.title('Variation within HEALpix, median=%4.2f mmag'%np.median(hp_std*1000))
        pyl.savefig('HPstd_hist.png',type='png')

    

    return starResults,HPResults


def printdict(dict):
    for key in dict.keys():
        if type(dict[key]) is int:
            print str(dict[key])+'    '+ key
        else:
            print '%4.3f'%(dict[key])+'    '+ key


if __name__ == "__main__":


    patchResults = patchPlots()


    
    starResults, HPResults = starPlots()

    print 'Patch Stats:'
    print '------------'
    printdict(patchResults)
    print '------------'

    print 'Star Stats:'
    print '------------'
    printdict(starResults)
    print '------------'

    print 'HEALpix binned results:'
    print '------------'
    printdict(HPResults)
    print '------------'
    


#To do

#maybe next do a set of properties (color distribution, stellar density, yadda yadda)
#come up with consistent xlim way of doing things
