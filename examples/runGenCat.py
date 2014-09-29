import numpy as np
from lsst.sims.selfcal.generation import genCatalog, offsets, visitOffsets
import lsst.sims.maf.db as db


# Read in an Opsim database
opsimDB = db.OpsimDatabase('sqlite:///opsimblitz2_1060_sqlite.db')
ralim=np.array([0,20])*np.pi/180.
declim=np.array([-5,-20])*np.pi/180.
cols = ['ditheredRA', 'ditheredDec', 'rotSkyPos', 'night', 'expMJD','fiveSigmaDepth' ]
visits = opsimDB.fetchMetricData(cols,
         '(fieldRA between %f and %f) and (fielddec between %f and %f)'%(ralim[0],ralim[1],declim[0],declim[1] ))
# Make dtype names more generic and add some other calcs:
visits = visitOffsets(visits, zpMax=1.)



