import os
import re
import lsst.sims.selfcal.analysis.color_correct as cc

cwd=os.getcwd()

name = re.sub('/.*/','',cwd)
tag='_'+name

cc.plotperstar(name, file='stats_rmsperstar.npz')
cc.bfpbsplot(name,bfnum='', from_main=True, file='stats_mobsmzpmbf.npz')
cc.checkColor(bestfit_star='test_bestfit_Star.dat', file='stats_mtruemmbf.npz')
cc.allStats(name,'stats_dmag')
cc.residCheck(name, fname='Mapall_dmag.png')
cc.residCheck(name, dmag='dmag_sin',fname='Mapsin_dmag.png')
cc.residCheck(name, dmag='dmag_flat_sin',fname='Mapflat_dmag.png')
#cc.residCheck(name,dmag='dmag_cloud_image+dmag_zp',visitID=1,fname='Mapcloudimage_dmag.png')
cc.residCheck(name,dmag='dmag_cloud_image',visitID=1,fname='Mapcloudimage_dmag.png')
cc.residCheck(name,dmag='dmag_cloud_image',fname='Mapcloudimage_night.png')
cc.nstarobs(name,file='stats_nobs.npz')
cc.rmsbins(name)
#cc.rmsrad(name)
