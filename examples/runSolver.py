from __future__ import print_function
from lsst.sims.selfcal.solver import lsqrSolver
from lsst.sims.selfcal.utils import fastRead
import numpy as np
import matplotlib.pylab as plt
from lsst.sims.selfcal.analysis.healplots import healbin
import healpy as hp



def robustRMS(val):
    iqr = np.percentile(val,75)-np.percentile(val,25)
    return iqr/1.349



mySolver = lsqrSolver(infile='observations.dat', patchOut='solvedPatch.npz', starOut='solvedStar.npz')

mySolver.run()


# now to compare the results!

starsFit = np.load('solvedStar.npz')['result']
trueStars = fastRead('starInfo.dat', dtype=zip(['starID','TrueMag', 'ra','dec'],
                                               [int,float,float,float]),
                                               delimiter=',')
trueStars.sort(order='starID')

# Remove any stars that were observed but not fit
mask = np.in1d(trueStars['starID'], starsFit['starID'], assume_unique=True)
trueStars = trueStars[mask]

resid = starsFit['fitMag'] - trueStars['TrueMag']

resid = resid-np.median(resid)

print('median fitMag - TrueMag = %f'%np.median(resid))
print('std (fitMag - TrueMag) = %f'%np.std(resid))
print('robust RMS = %f'%robustRMS(resid))
rrms = robustRMS(resid)

plt.hist(resid, bins=100, range=(-4*rrms,4*rrms))
plt.xlabel('Fit-True (mags)')
plt.ylabel('#')

plt.savefig('resid_hist.png')

healmap = healbin(np.radians(trueStars['ra']), np.radians(trueStars['dec']), resid)
hp.mollview(healmap, min =-4*rrms, max=4*rrms, unit='Residuals (mag)')
plt.savefig('resid_map.png')

