import numpy as np
from lsst.sims.selfcal.generation import genCatalog, offsets, visitOffsets
import lsst.sims.maf.db as db


lsstFilter = 'r'

# Read in an Opsim database
opsimDB = db.OpsimDatabase('sqlite:///opsimblitz2_1060_sqlite.db')
ralim=np.array([0,20])*np.pi/180.
declim=np.array([0,-20])*np.pi/180.
cols = ['ditheredRA', 'ditheredDec', 'rotSkyPos', 'night', 'expMJD','fiveSigmaDepth' ]

# Not sure why this query sucks...
#visits = opsimDB.fetchMetricData(cols,'(ditheredRA between %f and %f) and (ditheredDec between %f and %f) and filter="%s"'%(
#        ralim[0],ralim[1],declim[0],declim[1],lsstFilter ))
visits = opsimDB.fetchMetricData(cols,'ditheredRA < %f and ditheredRA > %f and ditheredDec > %f and ditheredDec < %f and  filter="%s"'%(ralim[1],ralim[0],declim[1],declim[0],lsstFilter  ))

# Make dtype names more generic and add any other stuff we want:
visits = visitOffsets.visitOffsets(visits, zpOff=1.)

offsetList=[]
# Systematic error floor for photometry
offsetList.append(offsets.OffsetSys() )
# SNR
offsetList.append(offsets.OffsetSNR() )

# Generate the catalog
genCatalog(visits, 'sqlite:///msrgb_1e6.sqlite', offsets=offsetList)

