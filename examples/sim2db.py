import os
import re
import lsst.sims.selfcal.analysis.color_correct as cc

cwd=os.getcwd()

name = re.sub('/.*/','',cwd)
tag='_'+name


cc.sim2db(name,tag)
