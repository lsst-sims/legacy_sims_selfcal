import numpy as np
from lsst.sims.selfcal.generation import genCatalog, offsets, visitOffsets
import lsst.sims.maf.db as db


lsstFilter = 'r'

# Read in an Opsim database
opsimDB = db.OpsimDatabase('sqlite:///opsimblitz2_1060_sqlite.db') #db.OpsimDatabase('sqlite:///opsimblitz2_1060_sqlite.db')
#ralim=np.array([0,360])*np.pi/180.
#declim=np.array([0,-90])*np.pi/180.
ralim=np.array([0,20])*np.pi/180.
declim=np.array([0,-20])*np.pi/180.
nightMax = 730

cols = ['ditheredRA', 'ditheredDec', 'rotSkyPos', 'night', 'expMJD','fiveSigmaDepth','obsHistID' ]

visits = opsimDB.fetchMetricData(cols,'ditheredRA < %f and ditheredRA > %f and ditheredDec > %f and ditheredDec < %f and  filter="%s" and night < %i'%(ralim[1],ralim[0],declim[1],declim[0],lsstFilter,nightMax  ))

# Make dtype names more generic and add any other stuff we want:
visits = visitOffsets(visits, zpOff=1.)

offsetList=[]
# Systematic error floor for photometry
offsetList.append(offsets.OffsetSys() )
# SNR
offsetList.append(offsets.OffsetSNR() )

# Generate the catalog
nobs, nstars = genCatalog(visits, 'sqlite:///msrgb_1e6.sqlite', offsets=offsetList,lsstFilter=lsstFilter)

