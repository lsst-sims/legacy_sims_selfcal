from lsst.sims.selfcal.solver import lsqrSolver
from lsst.sims.selfcal.utils import fastRead
import numpy as np
import matplotlib.pylab as plt

mySolver = lsqrSolver(infile='observations.dat', patchOut='solvedPatch.npz', starOut='solvedStar.npz')

mySolver.run()


# now to compare the results!

starsFit = np.load('solvedStar.npz')['result']
trueStars = fastRead('starInfo.dat', dtype=zip(['starID','TrueMag'], [int,float]),
                                               delimiter=',')
trueStars.sort(order='starID')

# Remove any stars that were observed but not fit
mask = np.in1d(trueStars['starID'], starsFit['starID'], assume_unique=True)
trueStars = trueStars[mask]

resid = starsFit['fitMag'] - trueStars['TrueMag']

resid = resid-np.median(resid)

print 'median fitMag - TrueMag = %f'%np.median(resid)
print 'std (fitMag - TrueMag) = %f'%np.std(resid)

import pdb ; pdb.set_trace()
