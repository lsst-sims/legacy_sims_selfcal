import numpy as np
from lsst.sims.selfcal.generation import genCatalog, offsets, visitOffsets
import lsst.sims.maf.db as db


lsstFilter = 'r'

# Read in an Opsim database
opsimDB = db.OpsimDatabase('sqlite:///opsimblitz2_1060_sqlite.db')
#ralim=np.array([0,360])*np.pi/180.
#declim=np.array([0,-90])*np.pi/180.
ralim=np.radians(np.array([40,60]) )
declim=np.radians(np.array([-40,-60]))
nightMax = 730

cols = ['ditheredRA', 'ditheredDec', 'rotSkyPos', 'night', 'expMJD','fiveSigmaDepth','obsHistID','transparency' ]

visits = opsimDB.fetchMetricData(cols,'ditheredRA < %f and ditheredRA > %f and ditheredDec > %f and ditheredDec < %f and  filter="%s" and night < %i'%(ralim[1],ralim[0],declim[1],declim[0],lsstFilter,nightMax  ))

# Make dtype names more generic and add any other stuff we want:
visits = visitOffsets(visits, zpOff=1.)


offsetList=[]
# Systematic error floor for photometry
offsetList.append(offsets.OffsetSys() )
# SNR
offsetList.append(offsets.OffsetSNR() )
# Clouds
offsetList.append(offsets.OffsetClouds() )

# Generate the catalog
nobs, nstars = genCatalog(visits, 'sqlite:///msrgb_1e6.sqlite', offsets=offsetList,lsstFilter=lsstFilter)


# Some simple timing tests.  Without Clouds turned on:
# time python runGenCat.py
# 7.989u 1.216s 1:32.52 9.9%      0+0k 1775128+40880io 37pf+0w
# time python runSolver.py
#LSQR finished
#The least-squares solution is good enough, given atol
#
#istop =       2   r1norm = 1.2e+02   anorm = 2.4e+04   arnorm = 2.6e-02
#itn   =     425   r2norm = 1.2e+02   acond = 2.1e+03   xnorm  = 1.5e+03
# std (fitMag - TrueMag) = 0.117947
#robust RMS = 0.000326
#8.848u 0.606s 0:17.25 54.7%     0+0k 41736+1968io 0pf+0w

# with clouds back on:
# time python runGenCat.py
# 708.105u 415.423s 20:02.94 93.3%        0+0k 1752512+43880io 0pf+0w
#time python runSolver.py
#median fitMag - TrueMag = 0.000000
#std (fitMag - TrueMag) = 0.117974
#robust RMS = 0.006536
#8.334u 0.527s 0:14.15 62.5%     0+0k 74248+1976io 72pf+0w
