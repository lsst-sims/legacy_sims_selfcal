import os
import re
import lsst.sims.selfcal.analysis.color_correct as cc

cwd=os.getcwd()

name = re.sub('/.*/','',cwd)
tag='_'+name

cc.flat_correct_nightblock(name, nx=40,ny=40)
#cc.flat_correct_nightblock(name, nx=80,ny=80)
