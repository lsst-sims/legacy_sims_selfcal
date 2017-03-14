from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import range
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import lsst.sims.selfcal.analysis.useful_input as ui
import pylab as pyl
import healpy as hp
import subprocess as sb
import multiprocessing
import os 
import lsst.sims.selfcal.analysis.catalog_plots_utils as cpu
from scipy import stats


def dictsort(indict,key):
    """for a dictionary of numpy arrays.  Normally, one would just use record arrays and np.sort(array,order=key) """
    order = np.argsort(indict[key])
    for key in list(indict.keys()):
        indict[key] = indict[key][order]
    return indict

def iqr(a):
    """interquartile range. 
    Should be even more robust to outliers than MAD."""
    s1 = stats.scoreatpercentile(a,25)
    s2 = stats.scoreatpercentile(a,75)
    return s2-s1

def MAD(a, axis=0):
    """ Median Absolute Deviation along given axis."""
    return np.median(np.fabs(a-np.median(a, axis=axis)), axis=axis)

def robustSigma(a, type='iqr'):
    """Compute a standard deviation that is robust to outliers.
    Choose between type='iqr' (interquartile range),
    and 'MAD' (Median Absolute Deviation)"""
    if type == 'iqr':
        c = 1.349
        sigma = iqr(a)/c
    if type == "MAD":
        c = 0.6745
        sigma = MAD(a)/c
    return sigma

def galac_b(ra,dec):
    """Compute galactic latitude """
    result = np.arcsin( np.sin(dec)*np.cos(62.87*np.pi/180.)-np.cos(dec)*np.sin(ra-282.86*np.pi/180.)*np.sin(62.87*np.pi/180))
    return result

def galac_l(ra,dec):
    #XXX there might be a typo in here...
    lcp = 123.932*np.pi/180
    alpha_GP = 192.85948*np.pi/180
    delta_GP = 27.12825*np.pi/180
    b = galac_b(ra,dec)
    lsin = np.cos(dec)*np.sin(ra-alpha_GP)/np.cos(b)
    lcos = np.cos(delta_GP)*np.sin(dec)-np.sin(delta_GP)*np.cos(dec)*np.cos(ra-alpha_GP)/np.cos(b)
    l = lcp-np.arctan2(lsin,lcos)
    if np.min(l) < 0:
        l[np.where(l < 0)] = l[np.where(l < 0)]+2.*np.pi
    return l


def print_obsfile(ofile, starobs_output="star_obs.dat", stars=None, subpatch=False):
    # This is the file that goes to selfCalib
    if ofile==None:
        ofile = open(starobs_output, 'w')
        #if subpatch:
        print("%s %s %s %s %s" %("#PatchID", "StarID", "StarObsMag", "StarMagObsErr", "FpPatchID"), file=ofile)
        #else:
        #    print >>ofile, "%s %s %s %s" %("#PatchID", "StarID", "StarObsMag", "StarMagObsErr")
    if stars!=None:
        for obj in range(0, len(stars['id'])):
            if subpatch:
                print("%d  %d  %f  %f %d" %(stars['fullpatch'][obj], stars['id'][obj],
                                                  stars['rmagobs'][obj], stars['magerr'][obj], stars['illum_patch'][obj]), file=ofile) #stars['subpatch'][obj])#
            else:
                print("%d  %d  %f  %f %d" %(stars['fullpatch'][obj], stars['id'][obj],
                                                  stars['rmagobs'][obj], stars['magerr'][obj], 0), file=ofile)
    return ofile

def split_obs(nside,shift_size=5., splitNS=False, NorthOnly=False, SouthOnly=False, blocksize=1e6):
    """ XXX kill this function XXX  replaced by just including neighboring healpixels"""
    """nside gives the number of sides for healpixels.  shift_size is in degrees and sets how  """


    #XXX-prob should put in a check so I don't sort repeatedly
    #if not os.path.isfile('star_obs.dat.s'):
    waitpop('sort -nk 1 star_obs.dat > star_obs.dat.s')


    npix = hp.nside2npix(nside)
    keys = ('visitid', 'patchid', 'subpatch', 'nstars', 'ra', 'dec') 
    patchfile='patchdata.dat'
    patchdat = ui.readDatafile(patchfile, keys)

    ord = np.argsort(patchdat['patchid'])
    for key in list(patchdat.keys()):  patchdat[key]=patchdat[key][ord]

    shift_size = shift_size*np.pi/180. #1.*np.pi/180.

    #compute galactic coords for healpixels and patches
    heal={}

    heal['ra'], heal['dec'] = hp.pix2ang(nside, np.arange(npix))
    heal['dec'] = heal['dec']-90.*np.pi/180.
    heal['b'] = galac_b(heal['ra'], heal['dec']) #np.arcsin(np.sin(heal['dec'])*np.cos(62.6*np.pi/180.)-np.cos(heal['dec'])*np.sin(heal['ra']-282.25*np.pi/180.)*np.sin(62.6*np.pi/180.))

    patchdat['dec'] = patchdat['dec']*np.pi/180.
    patchdat['ra'] = patchdat['ra']*np.pi/180.

    patchdat['b'] =  galac_b(patchdat['ra'], patchdat['dec']) 
    patchdat['l'] = galac_l(patchdat['ra'], patchdat['dec'])

    patchdat['dec']=patchdat['dec']+90.*np.pi/180. #stupid healpix wants 0-180 I think...


    patch_healpix = hp.ang2pix(nside, patchdat['dec'], patchdat['ra'])

    pyl.figure()
    cpu.plot_density(patchdat['ra'], patchdat['dec'], z=patch_healpix)
    pyl.savefig('heal_check.png', type='png')

    
    patch_healpix_shifted1 = hp.ang2pix(nside, patchdat['dec']+shift_size, patchdat['ra']+shift_size*2)
    patch_healpix_shifted2 = hp.ang2pix(nside, patchdat['dec']-shift_size, patchdat['ra']+shift_size*2)
    patch_healpix_shifted3 = hp.ang2pix(nside, patchdat['dec']+shift_size, patchdat['ra']-shift_size*2)
    patch_healpix_shifted4 = hp.ang2pix(nside, patchdat['dec']-shift_size, patchdat['ra']-shift_size*2)

    pyl.figure()
    cpu.plot_density(patchdat['ra'], patchdat['dec'], z=patch_healpix_shifted1)
    pyl.savefig('heal_check1.png', type='png')


    #ok, just go through and set duplicates to -1.  Also set any of the shifted ones that are not in patch_healpix to -1
    uhp = np.unique(patch_healpix)
    patch_healpix_shifted1[np.invert(np.in1d(patch_healpix_shifted1, uhp))] = -1
    patch_healpix_shifted2[np.invert(np.in1d(patch_healpix_shifted2, uhp))] = -1
    patch_healpix_shifted3[np.invert(np.in1d(patch_healpix_shifted3, uhp))] = -1
    patch_healpix_shifted4[np.invert(np.in1d(patch_healpix_shifted4, uhp))] = -1

    diff = patch_healpix_shifted1-patch_healpix
    good = np.where(diff == 0)
    patch_healpix_shifted1[good] = -1
    diff1 = patch_healpix_shifted2-patch_healpix
    diff2 = patch_healpix_shifted2-patch_healpix_shifted1
    good = np.where( (diff1 == 0) | (diff2 == 0))
    patch_healpix_shifted2[good] = -1
    diff1 = patch_healpix_shifted3-patch_healpix
    diff2 = patch_healpix_shifted3-patch_healpix_shifted1
    diff3 = patch_healpix_shifted3-patch_healpix_shifted2
    good = np.where( (diff1 == 0) | (diff2 == 0) | (diff3 == 0))
    patch_healpix_shifted3[good] = -1
    diff1 = patch_healpix_shifted4-patch_healpix
    diff2 = patch_healpix_shifted4-patch_healpix_shifted1
    diff3 = patch_healpix_shifted4-patch_healpix_shifted2
    diff4 = patch_healpix_shifted4-patch_healpix_shifted3
    good = np.where( (diff1 == 0) | (diff2 == 0) | (diff3 == 0) | (diff4 == 0))
    patch_healpix_shifted4[good] = -1
        
    pyl.figure()
    cpu.plot_density(patchdat['ra'], patchdat['dec'], z=patch_healpix_shifted1)
    pyl.savefig('heal_check1_post.png', type='png')

    heal_pixels = np.unique(patch_healpix)
    keys = [patch_healpix, patch_healpix_shifted1, patch_healpix_shifted2,patch_healpix_shifted3,patch_healpix_shifted4]
    filenames = ['hp1.dat','hp2.dat','hp3.dat','hp4.dat','hp5.dat']

    #ok, if a healpix stradles b =0, I want to split it in two
    # I could also just loop through each HP.  if the HP has patches w/ b> 0 and b< 0 and min l < xxx or max l < xxx, then break it apart.
    if splitNS:
        break_up = np.where( (patchdat['b'] < 0) & ( (patchdat['l'] < 120.*np.pi/180.) | (patchdat['l'] > 240)*np.pi/180.) & (np.abs(patchdat['b']) < 5*np.pi/180.))
        for key in keys:
            mask=np.where(key == -1)
            key[break_up] = key[break_up]+max(heal_pixels)+1
            key[mask] = -1
    
        


    print('reading star_obs.dat.s')
    star_obs = open('star_obs.dat.s', 'r')
    #read in just the patch id's from star_obs
    all_patch_ids = []
    for line in star_obs:
        if line.startswith('#'):
            continue
        linevalues = line.split()
        all_patch_ids.append(linevalues[0])

    all_patch_ids = np.array(all_patch_ids, dtype=np.int)
    star_obs.close()

    for j in np.arange(4):
        filename = filenames[j]
        patch_healpix = keys[j]
        hpid = all_patch_ids*0
    
        left = np.searchsorted(all_patch_ids, patchdat['patchid'])
        right = np.searchsorted(all_patch_ids, patchdat['patchid'], side='right')
        for i in np.arange(np.size(left)):
            hpid[left[i]:right[i]] =  patch_healpix[i]

        print('left %i, hpid %i, all_patch_ids %i, patch_healpix %i '%(np.size(left), np.size(hpid), np.size(all_patch_ids), np.size( patch_healpix)))
        file = open(filename,'w')
        print("#", file=file)
        for value in hpid:  print(value, file=file)
        file.close()
        #print j, np.size(hpid), np.size(all_patch_ids)
    
    waitpop('paste hp1.dat star_obs.dat.s > hp1_star_obs.dat')
    waitpop('paste hp2.dat star_obs.dat.s > hp2_star_obs.dat')
    waitpop("sed -i '/^-1/d' hp2_star_obs.dat" )
    waitpop('paste hp3.dat star_obs.dat.s > hp3_star_obs.dat')
    waitpop("sed -i '/^-1/d' hp3_star_obs.dat" )
    waitpop('paste hp4.dat star_obs.dat.s > hp4_star_obs.dat')
    waitpop("sed -i '/^-1/d' hp4_star_obs.dat" )
    waitpop('paste hp5.dat star_obs.dat.s > hp5_star_obs.dat')
    waitpop("sed -i '/^-1/d' hp5_star_obs.dat" )
    #waitpop("mv temp.dat hp2_star_obs.dat")
    waitpop('cat hp2_star_obs.dat >> hp1_star_obs.dat')
    waitpop('cat hp3_star_obs.dat >> hp1_star_obs.dat')
    waitpop('cat hp4_star_obs.dat >> hp1_star_obs.dat')
    waitpop('cat hp5_star_obs.dat >> hp1_star_obs.dat')
    waitpop(''' awk '{print $2" "$3" "$4" "$5 " "$6> ("h"$1"star_obs.dat")}' hp1_star_obs.dat''')
    waitpop('rm h#star_obs.dat')
    waitpop('mv hp1_star_obs.dat star_obs_hp1.dat')
    waitpop('mv hp2_star_obs.dat star_obs_hp2.dat')
    waitpop('mv hp3_star_obs.dat star_obs_hp3.dat')
    waitpop('mv hp4_star_obs.dat star_obs_hp4.dat')
    waitpop('mv hp5_star_obs.dat star_obs_hp5.dat')

def waitpop(command):
    print(command)
    result = sb.Popen(command, shell=True).wait()
    return result

def run_solver(solver, cpu_count=None):
    """Takes a string that is the path to the solver to run """
    files = sb.Popen('ls -s h*star_obs.dat', stdout=sb.PIPE, shell=True).stdout.read().split()
    
    n=len(files)/2
    sizes = files[0::2]
    sizes = np.array(sizes,dtype=np.float)
    files = np.array(files[1::2])
    ord = np.argsort(sizes)[::-1]
    files = files[ord] #run in order of largest to smallest as crude load-balancing

    if cpu_count == None:
        npool = multiprocessing.cpu_count()
    else:
        npool = cpu_count
    if (cpu_count != None) & (cpu_count < 0):
        npool = multiprocessing.cpu_count()+cpu_count

    pool = multiprocessing.Pool(processes=npool)

    commands=[]
    for i in np.arange(n):  commands.append('%s %s %s 1e-12 40000 > /dev/null'%(solver,files[i],files[i][0:-12]))

    results = pool.map_async(waitpop, commands,1).wait()
    pool.close()


def patch_combine():
    files = sb.Popen('ls h*bestfit_Patch.dat', stdout=sb.PIPE, shell=True).stdout.read().split()
    npatch=np.size(files)
    patchdat={}
    final_keys = ('patchid', 'patchmag', 'healid')
    for key in final_keys:    
        patchdat[key]=[]

    keys=('patchid','patchmag')

    for i in np.arange(npatch):
        temp = ui.readDatafile(files[i], keys, keytypes=('u8','float'))
        heal_num = np.float(files[i][1:-18])
        patchdat['patchid'].append(np.ravel(temp['patchid']))
        patchdat['healid'].append(np.ravel(temp['patchid']*0+heal_num))
        patchdat['patchmag'].append(np.ravel(temp['patchmag']))

    for key in list(patchdat.keys()):
        patchdat[key] = np.concatenate(patchdat[key])

    ord = np.argsort(patchdat['healid'])
    for key in list(patchdat.keys()):
        patchdat[key] = patchdat[key][ord]

    file = open('patch_obs.dat', 'w')
    print("%s %s %s %s %s" %("#HealID", "PatchID", "PatchMag", "PatchMagErr", "FpPatchID"), file=file)

    for obj in range(0,len(patchdat['patchid'])):
        print("%d  %d  %f  %f %d"%(patchdat['healid'][obj],patchdat['patchid'][obj], patchdat['patchmag'][obj], 0.001, 0), file=file)

    file.close()

def finalPatch(solver):
    files = sb.Popen('ls h*bestfit_Patch.dat', stdout=sb.PIPE, shell=True).stdout.read().split()
    #run the solver
    sb.Popen('%s patch_obs.dat patch 1e-12 40000 > /dev/null'%solver, shell=True).wait()
    sb.Popen('cp patch_bestfit_FpPatch.dat test_bestfit_FpPatch.dat', shell=True).wait()

    npatch = np.size(files)
    patchdat={}
    final_keys = ('patchid', 'patchmag', 'healid')
    for key in final_keys:    
        patchdat[key]=[]

    keys=('patchid','patchmag')

    
    for i in np.arange(npatch):
        temp = ui.readDatafile(files[i], keys, keytypes=('u8','float'))
        heal_num = np.float(files[i][1:-18])
        #insert if statement to deal with stupid patches that are not continuous
        patchdat['patchid'].append(np.ravel(temp['patchid']))
        patchdat['healid'].append(np.ravel(temp['patchid']*0+heal_num))
        patchdat['patchmag'].append(np.ravel(temp['patchmag']))

    for key in list(patchdat.keys()):
        patchdat[key] = np.concatenate(patchdat[key])

    ord = np.argsort(patchdat['healid'])
    for key in list(patchdat.keys()):
        patchdat[key] = patchdat[key][ord]

    keys=('hpid','hpzp')
    healdat = ui.readDatafile('patch_bestfit_Patch.dat',keys, keytypes=('u8','float'))
    keys=('patchid','patchmag')
    fit_patchdat = ui.readDatafile('patch_bestfit_Star.dat',keys, keytypes=('u8','float'))

    left = np.searchsorted(patchdat['healid'],healdat['hpid'])
    right = np.searchsorted(patchdat['healid'],healdat['hpid'], side='right')

    for i in np.arange(np.size(right)):
        patchdat['patchmag'][left[i]:right[i]] =  patchdat['patchmag'][left[i]:right[i]]-healdat['hpzp'][i]

    patchid = np.unique(patchdat['patchid'])
    patchmag = patchid*0
    ord = np.argsort(patchdat['patchid'])
    patchdat['patchmag'] = patchdat['patchmag'][ord]
    patchdat['patchid'] = patchdat['patchid'][ord]

    left = np.searchsorted(patchdat['patchid'],patchid)
    patchmag = patchdat['patchmag'][left]

    #now plug in the best fit mags where they exist
    left = np.searchsorted(patchid, fit_patchdat['patchid'])
    patchmag[left]=fit_patchdat['patchmag']

    file = open('test_bestfit_Patch.dat','w')

    for i in np.arange(np.size(patchmag)):
        print("%d %f"%(patchid[i], patchmag[i]), file=file)

    file.close()


def illum_combine():
    """combine the patch AND illumination corrections """
    files = sb.Popen('ls h*bestfit_Patch.dat', stdout=sb.PIPE, shell=True).stdout.read().split()
    #illums = sb.Popen('ls h*bestfit_FpPatch.dat', stdout=sb.PIPE, shell=True).stdout.read().split()
    #so_files = sb.Popen('ls h*bestfit_FpPatch.dat', stdout=sb.PIPE, shell=True).stdout.read().split()
    npatch=np.size(files)
    patchdat={}
    final_keys = ('patchid', 'patchmag', 'patchmagerr', 'illumid', 'healid')
    illum_keys = ('illumid','illummag')
    so_keys = ('patchid','starid','starmag','starmagerr','illumid','hpid')
    
    for key in final_keys:    
        patchdat[key]=[]

    keys=('patchid','patchmag')

    for i in np.arange(npatch):
        temp = ui.readDatafile(files[i], keys, keytypes=('u8','float'))
        heal_num = np.int(files[i][1:-18])
        temp_so = ui.readDatafile('h%u'%heal_num+'star_obs.dat', so_keys, keytypes=('u8','u8','float','float','u8','u8'))
        temp_illum = ui.readDatafile('h%u'%heal_num+'_bestfit_FpPatch.dat', illum_keys, keytypes=('u8','float'))
        #now I need to construct a line for each unique patchid+illumid pair
        temp_so['patch_plus_illum'] = np.core.defchararray.add(temp_so['patchid'].astype('|S10'),np.core.defchararray.add(np.core.defchararray.replace(np.zeros(np.size(temp_so['patchid']), dtype='|S1'), '',','), temp_so['illumid'].astype('|S10')))
        temp_so = dictsort(temp_so,'patch_plus_illum')
        upatch_plus_illum, unique_index = np.unique(temp_so['patch_plus_illum'],return_index=True )
        left = np.searchsorted(temp_so['patch_plus_illum'], upatch_plus_illum)
        right = np.searchsorted(temp_so['patch_plus_illum'], upatch_plus_illum, side='right')
        
        for key in so_keys:
            temp_so[key] = temp_so[key][unique_index]
        temp_so['nstars'] = right-left
        temp_so = np.core.records.fromarrays([temp_so['patchid'], temp_so['illumid'], temp_so['patchid']*0+heal_num, temp_so['starmag']*0,1./np.sqrt(temp_so['nstars']),temp_so['starmag']*0 ], names=('patchid','illumid','healid','patchmag', 'patchmagerr','illummag'))
        ord = np.argsort(temp['patchid'])
        for key in keys:
            temp[key] = temp[key][ord]
        temp_so.sort(order='patchid')
        left = np.searchsorted(temp_so['patchid'], temp['patchid'])
        right = np.searchsorted(temp_so['patchid'], temp['patchid'], side='right')
        for i in np.arange(np.size(left)):
            temp_so['patchmag'][left[i]:right[i]] = temp['patchmag'][i]

        temp_so.sort(order='illumid')
        ord = np.argsort(temp_illum['illumid'])
        for key in illum_keys:
            temp_illum[key] = temp_illum[key][ord]
        left = np.searchsorted(temp_so['illumid'], temp_illum['illumid'])
        right = np.searchsorted(temp_so['illumid'], temp_illum['illumid'], side='right')
        for i in np.arange(np.size(left)):
            temp_so['patchmag'][left[i]:right[i]] = temp_so['patchmag'][left[i]:right[i]] + temp_illum['illummag'][i]
        patchdat['patchid'].append(np.ravel(temp_so['patchid']))
        patchdat['healid'].append(np.ravel(temp_so['healid']))
        patchdat['patchmag'].append(np.ravel(temp_so['patchmag']))
        patchdat['illumid'].append(np.ravel(temp_so['illumid']))
        patchdat['patchmagerr'].append(np.ravel(temp_so['patchmagerr']))

    for key in list(patchdat.keys()):
        patchdat[key] = np.concatenate(patchdat[key])

    ord = np.argsort(patchdat['healid'])
    for key in list(patchdat.keys()):
        patchdat[key] = patchdat[key][ord]

    file = open('patch_obs.dat', 'w')
    print("%s %s %s %s %s" %("#HealID", "PatchID", "PatchMag", "PatchMagErr", "FpPatchID"), file=file)

    for obj in range(0,len(patchdat['patchid'])):
        print("%d  %d  %f  %f %d"%(patchdat['healid'][obj],patchdat['patchid'][obj], patchdat['patchmag'][obj], patchdat['patchmagerr'][obj], patchdat['illumid'][obj]), file=file)
        #print >>file, "%d  %d  %f  %f %d"%(patchdat['healid'][obj],patchdat['patchid'][obj], patchdat['patchmag'][obj], .01, patchdat['illumid'][obj])

    file.close()    
    
def finalIllum(solver):
    files = sb.Popen('ls h*bestfit_Patch.dat', stdout=sb.PIPE, shell=True).stdout.read().split()
    #run the solver
    sb.Popen('%s patch_obs.dat patch 1e-12 40000 > /dev/null'%solver, shell=True).wait()
    sb.Popen('cp patch_bestfit_Star.dat test_bestfit_Patch.dat', shell=True).wait()
    sb.Popen('cp patch_bestfit_FpPatch.dat test_bestfit_FpPatch.dat', shell=True).wait()

def finalStar_bad():
    """apply the healpix zeropoints to each block and collect all the bestfit star mags.
    I need to check that the multiple solutions per star are reasonable... is this an FpPatch issue?"""
    stardat={}
    keys = ('starid','starmag')
    for key in keys:
        stardat[key] = []
    healzp = ui.readDatafile('patch_bestfit_Patch.dat',('healid','healmag'), keytypes=('u8','float'))
    for i in np.arange(np.size(healzp['healid'])):
        temp = ui.readDatafile('h%i_bestfit_Star.dat'%healzp['healid'][i], keys, keytypes=('u8','float'))
        fp = ui.readDatafile('h%i_bestfit_FpPatch.dat'%healzp['healid'][i], ('id','mag'))
        temp['starmag'] = temp['starmag'] - healzp['healmag'][i] - fp['mag']
        stardat['starid'].append(np.ravel(temp['starid']))
        stardat['starmag'].append(np.ravel(temp['starmag']))
    for key in list(stardat.keys()):
        stardat[key] = np.concatenate(stardat[key])

    #Crap, looks like there's more jitter solution-to-solution than I expected.  might need to go back through the star_obs.dat file after all
    #convert stardat from dict of np arrays to a numpy rec array
    #sd = np.array(zip(stardat['starid'], stardat['starmag']), dtype=[('starid', 'u8'), ('starmag','float')] )
    #sd = np.core.records.fromarrays(stardat['starid'], stardat['starmag'],names='starid,starmag' )
    #sd.sort(order='starid')
    order = np.argsort(stardat['starid'])
    for key in list(stardat.keys()):
        stardat[key] = stardat[key][order]
    uid = np.unique(stardat['starid'])
    left = np.searchsorted(stardat['starid'],uid)
    right = np.searchsorted(stardat['starid'],uid, side='right')
    mags = np.zeros(left.size, dtype='float')
    stds = mags.copy()
    for i in np.arange(left.size):
        mags[i] = stardat['starmag'][left[i]:right[i]].mean()
        stds[i] = stardat['starmag'][left[i]:right[i]].std()
    file = open('test_bestfit_Star.dat', 'w')
    for i in np.arange(uid.size):
        print("%i %f"%(uid[i], mags[i]), file=file)
    file.close()
#crap, that didn't work well...looks like I have to go back to the star_obs file to construct the

def finalStar(infile='star_obs.dat', outfile='test_bestfit_Star.dat'):
    """read in the best fit patches, and the raw observations.  Apply patch zps to all the star magnitudes, then take a weighted mean of each star to make a bestfit star magnitude.  Also calc an stdev for each star"""
    patchzp = ui.readDatafile('test_bestfit_Patch.dat', ('patchid','patchzp'), keytypes=('u8','float') )
    star_obs = ui.readDatafile(infile, ('patchid','starid','starmag','starerr'),
                               keytypes=('u8','u8','float','float'))
    order = np.argsort(star_obs['patchid'])
    for key in list(star_obs.keys()):
        star_obs[key] = star_obs[key][order]
    left = np.searchsorted(star_obs['patchid'], patchzp['patchid'])
    right = np.searchsorted(star_obs['patchid'], patchzp['patchid'], side='right')
    for i in np.arange(left.size):
        star_obs['starmag'][left[i]:right[i]] = star_obs['starmag'][left[i]:right[i]] - patchzp['patchzp'][i]
    order = np.argsort(star_obs['starid'])
    for key in list(star_obs.keys()):
        star_obs[key] = star_obs[key][order]
    ustar = np.unique(star_obs['starid'])
    mag = np.zeros(ustar.size, dtype='float')
    magerr = np.zeros(ustar.size, dtype='float')
    magstd_raw = np.zeros(ustar.size, dtype='float')
    left = np.searchsorted(star_obs['starid'], ustar)
    right = np.searchsorted(star_obs['starid'], ustar, side='right')
    for i in np.arange(left.size):
        mag[i] = np.sum(star_obs['starmag'][left[i]:right[i]]/star_obs['starerr'][left[i]:right[i]]**2)/np.sum(1./star_obs['starerr'][left[i]:right[i]]**2)
        magerr[i] = 1./np.sqrt(np.sum(1./star_obs['starerr'][left[i]:right[i]]**2))
        magstd_raw[i] = star_obs['starmag'][left[i]:right[i]].std()
    file = open(outfile, 'w')
    print('#StarID   Mag    Err   STDEV', file=file)
    for i in np.arange(ustar.size):
        print("%i %f %f %f"%(ustar[i], mag[i], magerr[i], magstd_raw[i]), file=file)
    file.close()



    
def finalStarP(infiles='h*star_obs.dat', outfile='test_bestfit_Star.dat', outstats='stats_stars'):
    """ read in all the healpix files and apply patch zeropoints.  Compute statistics for each star as well as final magnitude """
    
    patchzp = ui.readDatafile('test_bestfit_Patch.dat', ('patchid','patchzp'), keytypes=('u8','float') )
    illumzp = ui.readDatafile('test_bestfit_FpPatch.dat', ('illumid', 'illumzp'), keytypes=('u8','float'))
    files = np.array(sb.Popen('ls %s'%infiles, stdout=sb.PIPE, shell=True).stdout.read().split())
    hps = files.copy()
    for i in np.arange(np.size(hps)):
        hps[i]=hps[i].replace('h','').replace('star_obs.dat','')
    hps=hps.astype('int')
    #true_stars = ui.readData('stardata.dat', ('starid','stardbid','ra','dec','mag','color'), keytypes=('u8','u8','float','float','float','float'))
    results = {}
    result_keys = ('starid','nobs','mag','stdev','stdev_IQR')
    for key in result_keys:
        results[key]=[]

    for i in np.arange(np.size(files)):
        #read in file
        star_obs = ui.readDatafile(files[i], ('patchid','starid','starmag','starerr','fppatch','hp'),
                                   keytypes=('u8','u8','float','float','u8','u8'))
        in_hp = np.where(star_obs['hp'] == hps[i])#crop down to just those stars in the center HP
        for key in list(star_obs.keys()):
            star_obs[key] = star_obs[key][in_hp]
        order = np.argsort(star_obs['patchid'])
        for key in list(star_obs.keys()):
            star_obs[key] = star_obs[key][order]

        left = np.searchsorted(star_obs['patchid'], patchzp['patchid'])
        right = np.searchsorted(star_obs['patchid'], patchzp['patchid'], side='right')
        for i in np.arange(left.size): #apply patch zeropoints
            star_obs['starmag'][left[i]:right[i]] = star_obs['starmag'][left[i]:right[i]] - patchzp['patchzp'][i]

        order = np.argsort(star_obs['fppatch'])
        for key in list(star_obs.keys()):
            star_obs[key] = star_obs[key][order]
        left = np.searchsorted(star_obs['fppatch'], illumzp['illumid'])
        right = np.searchsorted(star_obs['fppatch'], illumzp['illumid'], side='right')
        for i in np.arange(left.size):
            star_obs['starmag'][left[i]:right[i]] = star_obs['starmag'][left[i]:right[i]] - illumzp['illumzp'][i]
        
        
        order = np.argsort(star_obs['starid'])
        for key in list(star_obs.keys()):
            star_obs[key] = star_obs[key][order]
        ustar = np.unique(star_obs['starid'])
        left = np.searchsorted( star_obs['starid'], ustar)
        right = np.searchsorted(star_obs['starid'], ustar,side='right')
        temp_results = {}
        for key in result_keys:
            temp_results[key] = np.zeros(left.size)
        temp_results['starid'] = ustar.copy()
        for i in np.arange(left.size):
            temp_results['nobs'][i] = np.size(star_obs['starmag'][left[i]:right[i]])
            temp_results['mag'][i] = np.sum(star_obs['starmag'][left[i]:right[i]]/star_obs['starerr'][left[i]:right[i]]**2)/np.sum(1./star_obs['starerr'][left[i]:right[i]]**2)
            temp_results['stdev'][i] = np.std(star_obs['starmag'][left[i]:right[i]])
            temp_results['stdev_IQR'][i] = robustSigma(star_obs['starmag'][left[i]:right[i]])
        for key in result_keys:
            results[key].append(temp_results[key].copy())
    for key in result_keys:
        results[key] = np.concatenate(results[key])
    ord = np.argsort(results['starid'])
    for key in result_keys:
        results[key] = results[key][ord]
        #output test_bestfit_Star w/ID, mag and then the rest.
        #I guess I could put the print inside the loop if I was worried about memory
    file = open('test_bestfit_Star.dat','w')
    for i in np.arange(np.size(results['starid'])):
        print('%d %f %f %f %d'%(results['starid'][i], results['mag'][i],
                                        results['stdev'][i], results['stdev_IQR'][i],results['nobs'][i]), file=file)
    file.close()

