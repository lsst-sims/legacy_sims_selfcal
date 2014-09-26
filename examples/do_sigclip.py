import lsst.sims.selfcal.analysis.color_correct as cc
import os
import re
cwd=os.getcwd()

name = re.sub('/.*/','',cwd)
tag='_'+name

cc.star_sigma_clip(name)
