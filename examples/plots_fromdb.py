import os
import re
import lsst.sims.selfcal.analysis.color_correct as cc

cwd=os.getcwd()

name = re.sub('/.*/','',cwd)
tag='_'+name

cc.corrected2db(name,tag,bfnum='')
cc.bfpbsplot(name,bfnum='')
cc.checkColor(bestfit_patch='test_bestfit_Patch.dat',bestfit_star='test_bestfit_Star.dat')



