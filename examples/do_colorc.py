import os
import re
import lsst.sims.selfcal.analysis.color_correct as cc

cwd=os.getcwd()

name = re.sub('/.*/','',cwd)
tag='_'+name


cc.makeDeltamap(name)
cc.applyDeltamap(name)
cc.color_corrected2db(name,'')
#cc.bfpbsplot(name,bfnum='')
#cc.checkColor(bestfit_patch='test_bestfit_Patch.dat',bestfit_star='test_bestfit_Star.dat')
cc.plotDelta()



