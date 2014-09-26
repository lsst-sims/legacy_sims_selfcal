#####
#  Lynne Jones, ljones@astro.washington.edu
#  svn version info : $Id$
#
#
#####


# general python modules
import numpy as n
import numpy.ma as M
import pylab as pyl

# self-cal python modules
import lsst.sims.selfcal.analysis.useful_input as ui
import lsst.sims.selfcal.analysis.catalog_plots_utils as cpu

deg2rad = n.pi/180.0
rad2deg = 180.0/n.pi
sqrt2 = n.sqrt(2)


#-- star magnitudes (true and calibrated / best fit)

def read_true_stars(stardatfile='stardata.dat'):
    """Read a stardata.dat file from simSelfCalib.py"""
    keys = ('id', 'dbid', 'ra', 'dec', 'magtrue', 'color')
    stardat = ui.readDatafile(stardatfile, keys)
    return stardat

def read_bestfit_stars(starcalfile='test_bestfit_Star.dat'):
    """Read a bestfit_Star.dat file from simple.x"""
    keys = ('id', 'magcal')
    starcal = ui.readDatafile(starcalfile, keys)
    return starcal

def match_stars_calvstrue(stardat, starcal, sorted=True):
    """Match up true information to calibration information, return extended starcal.
    Note that stardat and starcal do not have to be the same length; in fact, starcal
    will usually be slightly shorter than stardat."""
    stardat['magcal'] = n.zeros(len(stardat['id']), dtype='float') - 99
    # starcal['id'] should be sorted (from simple.x)! 
    # stardat['id'] should be sorted (although stardat['dbid'] is probably not)!
    # stardat['id'] is the same identification number as starcal['id'] .. to track back
    # to the database, however, stardat['dbid'] would be necessary. 
    if sorted: 
        j = 0
        for i in range(len(starcal['id'])):
            # Step through (sorted?) stardat looking for stars which match the starcal ids.
            # stardat presumed sorted, but not necessarily
            # starcal may skip some stardat values
            while (starcal['id'][i] != stardat['id'][j]) & (j<len(stardat['id'])):
                j += 1
            if j> len(stardat['id']):
                sorted = False
            stardat['magcal'][j] = starcal['magcal'][i]
    if not sorted:
        for i in range(len(starcal['id'])):
            idx = n.where(stardat['id']==starcal['id'][i])
            idx = idx[0][0]
            stardat['magcal'][idx] = starcal['magcal'][i]
    condition = (stardat['magcal']!=-99)
    for key in stardat.keys():
        stardat[key] = stardat[key][condition]
    return stardat

def adjust_stars_calvstrue(stardat):
    """Adjust the calibrated magnitudes to the true magnitudes .. ie. remove the floating zp."""
    stardat['magdiff'] = stardat['magtrue'] - stardat['magcal']
    meanin = stardat['magdiff'].mean()
    subset = stardat['magdiff'][n.where((abs(stardat['magdiff']-meanin)<3) & (stardat['magcal']!=-99))]
    floatzp = subset.mean()
    stardat['magdiff'] = stardat['magdiff'] - floatzp
    meanout = stardat['magdiff'].mean()
    rms = stardat['magdiff'].std()
    print "Before adjustment: mean %f.  After adjustment: mean %f with stdev %f" %(meanin, meanout, rms)
    #return stardat, meanin, meanout, rms
    return stardat
    
def plot_stars_dmagcaltrue(stardat, maglim=None, etitle="", newfig=True,
                           starcal=None, z_method=n.mean):
    adjust_stardat = False
    if maglim !=0:
        maglim=n.round(n.array(maglim)*1000, decimals=1)
    if 'magcal' not in stardat.keys():
        if starcal == None:
            raise Exception("Must either pass in starcal as well, or match_calvstrue previously")
        stardat = match_calvstrue(stardat, starcal)
        adjust_stardat = True
    if 'magdiff' not in stardat.keys():
        stardat = adjust_calvstrue(stardat)
        adjust_stardat = True
    if 'x' not in stardat.keys():
        stardat['x'], stardat['y'] = cpu.hammer_project_toxy(stardat['ra']*deg2rad,
                                                            stardat['dec']*deg2rad)
        adjust_stardat = True
    gridsize = cpu.calc_gridsize(stardat, binsize='patch')
    cpu.plot_density(stardat['x'], stardat['y'], stardat['magdiff']*1000, z_method=z_method,
                     zlimits=maglim, gridsize=gridsize, radecgrid=True, newfig=newfig, cb_label='mmag')
    figtitle = "Stars dMag(true-bestfit) " + etitle
    pyl.title(figtitle, fontsize='medium')
    return stardat

def histogram_stars_dmagcaltrue(stardat, nbins=100, range=None, etitle="", label='', newfig=True):
    """Histogram the magnitude error."""
    if range != None:
        range=n.array(range)*1000
    if newfig:
        pyl.figure()
    if label!='':
        num, b, p = pyl.hist(stardat['magdiff']*1000, bins=nbins, range=range, label=label,
                           facecolor=None, histtype='step')
        pyl.legend(numpoints=1, loc='upper right')
    else:
        num, b, p = pyl.hist(stardat['magdiff']*1000, bins=nbins, range=range, histtype='bar')
    pyl.xlabel("dMag(true-bestfit) (mmag)")
    string = "Stars dMag(true-bestfit)" + etitle
    pyl.title(string, fontsize='medium')
    return

#--- visits

def read_visitfile(visitfile='visit.dat'):
    """Read the visit information from simSelfCalib.py"""
    keys = ('visitid', 'ra', 'dec', 'skyrot', 'time', 'zpOff', 'zpGrad', 'zpGradDir',
            'colorterm', 'colortermRad')
    visits = ui.readDatafile(visitfile, keys)
    return visits

def plot_visits(visits, visitlim=None, etitle="", newfig=True):
    """Plot the location of the visits across the sky. Color code by number of visits."""
    if 'x' not in visits.keys():
        visits['x'], visits['y'] = cpu.hammer_project_toxy(visits['ra']*deg2rad,
                                                          visits['dec']*deg2rad)
    gridsize = cpu.calc_gridsize(visits,binsize='fov', rad_fov=1.8)
    cpu.plot_density(visits['x'], visits['y'], gridsize=gridsize, zlimits=visitlim,
                    radecgrid=True, newfig=newfig)#, zscale='log')
    pyl.title("Visit density " + etitle, fontsize='medium')
    return

#---- patches

def read_true_patchfile(patchfile='patchdata.dat'):
    """Read the patch file information from simSelfCalib.py"""
    keys = ('visitid', 'patchid', 'subpatch', 'nstars', 'ra', 'dec', 'magtrue', 'color',
            'dmag', 'dmag_rms', 'dmag_rand', 'dmag_rand_rms', 'dmag_zp', 'dmag_zp_rms',
            'dmag_color', 'dmag_color_rms')
    patchdat = ui.readDatafile(patchfile, keys)
    # dmag = stars' (true magnitude - 'observed') magnitude, average 
    #patchdat['dmag'] = -1*patchdat['dmag']
    return patchdat

def read_bestfit_patchfile(bestfitpatchfile='test_bestfit_Patch.dat'):
    """Read the best fit patch file information from simple.x"""
    keys = ('patchid', 'dmagcal')
    patchcal = ui.readDatafile(bestfitpatchfile, keys)
    # dmagcal = zeropoint calculated by calsim
    return patchcal

def plot_nvisits_patch(patchdat, subpatch=None, patchlim=None, newfig=True):
    """Plot the location of the (sub)patches across the sky. Color code by number of visits.
    Also plot the number of stars observed in each patch. """
    if 'x' not in patchdat.keys():
        patchdat['x'], patchdat['y'] = cpu.hammer_project_toxy(patchdat['ra']*deg2rad,
                                                            patchdat['dec']*deg2rad)
    if subpatch != None:
        condition = (patchdat['subpatch'] == subpatch)
        x = patchdat['x'][condition]
        y = patchdat['y'][condition]
    else:
        x = patchdat['x']
        y = patchdat['y']
    gridsize = cpu.calc_gridsize(patchdat, binsize='patch')
    cpu.plot_density(x, y, gridsize=gridsize, zlimits=patchlim, radecgrid=True, newfig=newfig)
    pyl.title("Patch Density", fontsize='medium')
    return

def histogram_nstars_patch(patchdat, show_subpatch=False, newfig=True):
    """Histogram the number of stars per patch."""
    if newfig:
        pyl.figure()
    nbins = 100
    if show_subpatch:
        pmin = patchdat['subpatch'].min()
        pmax = patchdat['subpatch'].max()
        for ipatch in range(pmin, pmax, 1):
            condition = (patchdat['subpatch']==ipatch)
            counts = patchdat['nstars'][condition]
            num, b, p = pyl.hist(counts, bins=nbins, label="%d" %(ipatch), histtype='step', facecolor=None)
        pyl.legend(numpoints=1, loc='upper right')
    else:
        num, b, p = pyl.hist(patchdat['nstars'], bins=nbins, normed=1)
    pyl.xlabel("Number of stars per patch")
    return
     
def plot_nstars_patch(patchdat, z_method=n.mean, zlimits=None, newfig=True):
    """Plot the n.min/mean/std number of stars ever seen in a particular patch."""
    if 'x' not in patchdat.keys():
        patchdat['x'], patchdat['y'] = cpu.hammer_project_toxy(patchdat['ra']*deg2rad,
                                                            patchdat['dec']*deg2rad)
    gridsize = cpu.calc_gridsize(patchdat, binsize='patch')
    cpu.plot_density(patchdat['x'], patchdat['y'], z=patchdat['nstars'], z_method=z_method,
                    zlimits=zlimits, gridsize=gridsize, radecgrid=True, newfig=newfig)
    pyl.title("Patch Minimum Nstars")
    return

def plot_minstars_patch(patchdat, minstars=2, subpatch=None, newfig=True):
    """Scatter plot the x/y locations of patchdat where the number of stars in a patch
    is less than minstars."""
    if newfig:
        pyl.figure()
    condition = (patchdat['nstars']<minstars)
    if subpatch!=None:
        condition = condition & (patch['subpatch']==subpatch)
    tempx = patchdat['x'][condition]
    tempy = patchdat['y'][condition]
    print "%f%s (%d) of the total %d patchdat have less than %f stars" \
          %(float(len(tempx))/len(patchdat['nstars'])*100.0, 
            "%s", len(tempx), len(patchdat['nstars']), minstars)
    pyl.plot(tempx, tempy, 'k.')
    cpu.make_radec_grid(cpu.hammer_project_toxy)
    string = "Locations of patchdat with less than %d stars" %(minstars)
    if subpatch!=None:
        string = string + ", for subpatch %d" %(subpatch)
    pyl.title(string, fontsize='medium') 
    return

def match_patch_calvstrue(patchdat, patchcal, sorted=True):
    """Match up true information to calibration information, return extended patch data.
    Note that patchdat and patchcal do not have to be the same length; in fact, patchcal
    will usually be slightly shorter than patchdat."""
    patchdat['dmagcal'] = n.zeros(len(patchdat['patchid']), dtype='float') - 99.0
    # patchcal['id'] should be sorted (from simple.x)! 
    # patchdat['id'] should be sorted too from simSelfCalib.py
    if sorted: 
        j = 0
        for i in range(len(patchcal['patchid'])):
            # Step through (sorted?) patchdat looking for patches which match the patchcal ids.
            # some patches (uncalibrateable) might be skipped
            while (patchcal['patchid'][i] != patchdat['patchid'][j]) & (j<len(patchdat['patchid'])):
                j += 1
            if j> len(patchdat['patchid']):
                sorted = False
            patchdat['dmagcal'][j] = patchcal['dmagcal'][i]
    if not sorted:
        for i in range(len(patchcal['patchid'])):
            idx = n.where(patchdat['patchid']==patchcal['patchid'][i])
            idx = idx[0][0]
            patchdat['dmagcal'][idx] = patchcal['dmagcal'][i]
    condition = (patchdat['dmagcal']!=-99)
    for key in patchdat.keys():
        patchdat[key] = patchdat[key][condition]
    return patchdat

def adjust_patch_calvstrue(patchdat):
    """Adjust the calibrated zeropoint magnitudes to the true zeropoint magnitudes ..
    ie. remove the floating zp."""
    patchdat['magdiff'] = patchdat['dmag'] - patchdat['dmagcal']
    meanin = patchdat['magdiff'].mean()
    subset = patchdat['magdiff'][n.where(abs(patchdat['magdiff']-meanin)<3)]
    floatzp = subset.mean()
    patchdat['magdiff'] = patchdat['magdiff'] - floatzp
    meanout = patchdat['magdiff'].mean()
    rms = patchdat['magdiff'].std()
    print "Before adjustment: mean %f.  After adjustment of %f, find mean %f with stdev %f" \
          %(meanin, floatzp, meanout, rms)
    #return stardat, meanin, meanout, rms
    return patchdat
    
def plot_patch_dmagcaltrue(patchdat, maglim=None, etitle="", newfig=True,
                           patchcal=None, z_method=n.mean):
    adjust_patchdat = False
    if maglim != None:
        maglim=n.round(n.array(maglim)*1000, decimals=1)
    if 'dmagcal' not in patchdat.keys():
        if patchcal == None:
            raise Exception("Must either pass in patchcal as well, or match_patch_calvstrue previously")
        patchdat = match_calvstrue(patchdat, starcal)
        adjust_patchdat = True
    if 'magdiff' not in patchdat.keys():
        patchdat = adjust_calvstrue(patchdat)
        adjust_patchdat = True
    if 'x' not in patchdat.keys():
        patchdat['x'], patchdat['y'] = cpu.hammer_project_toxy(patchdat['ra']*deg2rad,
                                                            patchdat['dec']*deg2rad)
        adjust_patchdat = True
    gridsize = cpu.calc_gridsize(patchdat, binsize='patch')
    cpu.plot_density(patchdat['x'], patchdat['y'], patchdat['magdiff']*1000., z_method=z_method,
                     zlimits=maglim, gridsize=gridsize, radecgrid=True, newfig=newfig, cb_label='mmag')
    figtitle = "Patch dMag(true-bestfit) " + etitle
    pyl.title(figtitle, fontsize='medium')
    return patchdat

def plot_patch_nstars(patchdat, maglim=None, etitle="", newfig=True,
                           patchcal=None, z_method=n.mean):
    adjust_patchdat = False
    if maglim != None:
        maglim=n.array(maglim)
    if 'dmagcal' not in patchdat.keys():
        if patchcal == None:
            raise Exception("Must either pass in patchcal as well, or match_patch_calvstrue previously")
        patchdat = match_calvstrue(patchdat, starcal)
        adjust_patchdat = True
    if 'magdiff' not in patchdat.keys():
        patchdat = adjust_calvstrue(patchdat)
        adjust_patchdat = True
    if 'x' not in patchdat.keys():
        patchdat['x'], patchdat['y'] = cpu.hammer_project_toxy(patchdat['ra']*deg2rad,
                                                            patchdat['dec']*deg2rad)
        adjust_patchdat = True
    gridsize = cpu.calc_gridsize(patchdat, binsize='patch')
    cpu.plot_density(patchdat['x'], patchdat['y'], patchdat['nstars'], z_method=z_method,
                     zlimits=maglim, gridsize=gridsize, radecgrid=True, newfig=newfig, cb_label='')
    figtitle = "N stars per patch" + etitle
    pyl.title(figtitle, fontsize='medium')
    return patchdat


def plot_patch_zpoffset(patchdat, maglim=None, etitle="", newfig=True, z_method=n.mean):
    if 'x' not in patchdat.keys():
        patchdat['x'], patchdat['y'] = cpu.hammer_project_toxy(patchdat['ra']*deg2rad,
                                                               patchdat['dec']*deg2rad)
    gridsize = cpu.calc_gridsize(patchdat, binsize='patch')
    cpu.plot_density(patchdat['x'], patchdat['y'], patchdat['dmag'], z_method=z_method,
                     zlimits=maglim, gridsize=gridsize, radecgrid=True, newfig=newfig, cb_label='mag')
    figtitle = "Zeropoint Offsets " + etitle
    pyl.title(figtitle, fontsize='medium')
    return 

def histogram_patch_dmagcaltrue(patchdat, nbins=100, range=None, etitle="", label=None, newfig=True):
    """Histogram the magnitude error."""
    if range != None:
        range=n.array(range)*1000
    if newfig:
        pyl.figure()
    if label!=None:
        num, b, p = pyl.hist(patchdat['magdiff']*1000, bins=nbins, range=range, label=label,
                           facecolor=None, histtype='step')
        pyl.legend(numpoints=1, loc='upper right')
    else:
        num, b, p = pyl.hist(patchdat['magdiff']*1000, bins=nbins, range=range, histtype='bar')
    pyl.xlabel("dMag(true-bestfit) (mmag)")
    string = "Patch dMag(true-bestfit) " + etitle
    pyl.title(string, fontsize='medium')
    return

def histogram_patch_zpoffset(patchdat, nbins=100, range=None, etitle="", label=None, newfig=True):
    """Histogram the magnitude error."""
    if newfig:
        pyl.figure()
    if label!=None:
        num, b, p = pyl.hist(patchdat['dmag'], bins=nbins, range=range, label=label,
                           facecolor=None, histtype='step')
        pyl.legend(numpoints=1, loc='upper right')
    else:
        num, b, p = pyl.hist(patchdat['dmag'], bins=nbins, range=range, histtype='bar')
    pyl.xlabel("Zeropoint Offsets (Mags)")
    string = "Zeropoint Offsets " + etitle
    pyl.title(string, fontsize='medium')
    return

#--- star observations (these can be large)

def read_starsobs(starobsfile='star_obs.dat'):
    """Read the star_obs.dat file from simSelfCalib.py that is input to simple.x"""
    keys = ('patchid', 'starid', 'rmagobs', 'magerr')
    starobs = ui.readDatafile(starobsfile, keys)
    return starobs

def read_masterobsfile(masterfile='master_cal.dat', readkeys=None, count_nobs_star=True):
    """Read selected information from the master data file. Since this file can be large,
    generally a single piece of information will want to be extracted. Multiple columns can
    be read at once though, if desired. In addition, some calculations involving counting
    (such as number of observations of a particular star) can be done while reading the file."""
    f = open(masterfile, 'r')
    keys = ('visitid', 'fullpatch', 'subpatch', 'starid', 'dbid', 'rmagobs', 'magerr',
            'dmag_var', 'dmag_snr', 'dmag_zp', 'dmag_zp_dist', 'dmag_color', 'dmag_color_dist',
            'dmag_sin', 'dMag_flat_sin','dMag_cloud', 'dMag_cloud_image','dMag_rr',
            'rmagtrue', 'color', 'ra', 'dec', 'ra_f', 'dec_f', 'X', 'Y', 'Night')
    keytypes = ('int', 'int', 'int', 'int', 'int', 'float', 'float',
                'float', 'float', 'float', 'float', 'float', 'float',
                'float', 'float', 'float', 'float','float', 'float', 'float',
                'float', 'float','float','float', 'float', 'float', 'float')
    if readkeys != None:
        readcols = []
        for i in readkeys:
            readcols.append(keys.index(i))
        readcols = n.array(readcols, dtype='int')
    # Count number of observations of each star, store in dictionary.
    if count_nobs_star:
        count_nobs={}
    # Read desired values from master obs file, line by line.
    data = {}
    for readkey in readkeys:
        data[readkey] = []
    for line in f:
        values = line.split()
        if line.startswith("!"):
            continue
        if line.startswith("#"):
            continue
        for i in range(0, len(readkeys)):
            data[readkeys[i]].append(values[readcols[i]])
        if count_nobs_star:
            starid = values[keys.index('starid')]
            if starid in count_nobs.keys():                
                count_nobs[starid] += 1
            else:
                count_nobs[starid] = 1
    for i in range(0, len(readkeys)):
        data[readkeys[i]] = n.array(data[readkeys[i]], dtype=keytypes[readcols[i]])
    if count_nobs_star:
        return data, count_nobs
    else:
        return data

def translate_count_nobs(count_nobs):
    nobs_star = n.zeros(len(count_nobs.keys()), dtype='int')
    i = 0
    for starid in count_nobs:
        nobs_star[i] = count_nobs[starid]
        i += 1
    return nobs_star

def histogram_starnobs(nobs_star, newfig=True, nbins=100, xlim=None, ylim=None):
    """Histogram the number of observations per star. """
    if newfig:
        pyl.figure()
    # Check to see if user sent in array (already translated) or if nobs_star needs translating. 
    counts = nobs_star
    try:
        counts[0]
    except KeyError:
        counts = translate_count_nobs(nobs_star)
    # Do histogram.
    num, b, p = pyl.hist(counts, bins=nbins, histtype='step', facecolor=None)
    if xlim!=None:
        pyl.xlim(xlim)
    if ylim!=None:
        pyl.ylim(ylim)
    pyl.xlabel("Number of observations per star")
    return

