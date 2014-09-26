import numpy as np
import matplotlib
matplotlib.use('Agg')
import lsst.sims.selfcal.generation.split_obs as so

#just use awk to split the observations onto healpatches
#so.waitpop('''awk 'NR==1{next} {print $1" "$2" "$3" "$4" "$5 > ("h"$6"star_obs.dat")}' star_obs.dat ''')
#so.waitpop('''awk 'NR==1{next} {if ($7 != -1) print $1" "$2" "$3" "$4" "$5 >> ("h"$7"star_obs.dat")}' star_obs.dat ''')
#so.waitpop('''awk 'NR==1{next} {if ($8 != -1) print $1" "$2" "$3" "$4" "$5 >> ("h"$8"star_obs.dat")}' star_obs.dat ''')
#so.waitpop('''awk 'NR==1{next} {if ($9 != -1) print $1" "$2" "$3" "$4" "$5 >> ("h"$9"star_obs.dat")}' star_obs.dat ''')

#better way that goes through the file only once
so.waitpop('rm h*star_obs.dat h*_bestfit_*.dat h*_restart.dat') #clean up any old files still wandering around.
so.waitpop(''' awk 'NR==1{next} {print $1" "$2" "$3" "$4" "$5" "$6>> ("h"$6"star_obs.dat") ; if ($7 != -1 && $7 != $6) print $1" "$2" "$3" "$4" "$5" "$6 >> ("h"$7"star_obs.dat") ; if ($8 != -1 && $8 != $6) print $1" "$2" "$3" "$4" "$5" "$6 >> ("h"$8"star_obs.dat") ; if ($9 != -1 && $9 != $6) print $1" "$2" "$3" "$4" "$5" "$6 >> ("h"$9"star_obs.dat") ; if ($10 != -1 && $10 != $6) print $1" "$2" "$3" "$4" "$5" "$6 >> ("h"$10"star_obs.dat")}' star_obs.dat ''' )

#XXX -- make the star_ids something reasonable in each healpix

#run the single-core solver on each healpixel.  
so.run_solver('$solver2', cpu_count = -1)

#combine the healpixel solutions
#so.patch_combine()
so.illum_combine()


#run the multi-core solver to tie the healpixels to the same floating zeropoint
#so.finalPatch('$solver')
so.finalIllum('$solver')


#loop back through the stars and apply the patches
so.finalStarP()

#plot up how the final patch solution look
#execfile('/astro/store/shared-scratch1/yoachim/lsst_devel_winter12/Linux64/sims/selfcal/analysis/trunk/examples/patch_simple.py')

#rm h*star_obs* h*_bestfit* h*_restart*  *.png
#python healsolve.py ; python $sce/patch_simple.py  
