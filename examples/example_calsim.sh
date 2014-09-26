#!/bin/csh

###Example script to be used if there is a color gradient that needs to be corrected for.


#setenv pp <path to your sandx>/lsst_devel/Linux64/selfcal/python/selfcal/analyze
setenv pp1 ${SELFCAL_GENERATION_DIR}/python/lsst/sims/selfcal/generation
setenv pp2 ${SELFCAL_ANALYSIS_DIR}/examples/

#######################
#Example run of the calsim.  Assumes all the code is in your path

#generate the observed magnitudes for the observations
#python $pp/../gendat/simSelfCalib.py simSelfCalib.input
python $pp1/simSelfCalib.py ${SELFCAL_DATA}/simSelfCalib.input

#run the calibration fitter on the observed stellar mags
${SELFCAL_SOLVER_DIR}/bin/simple.x star_obs.dat test 1.e-7 10000

#make some plots
python $pp2/simple_plots_star.py
python $pp2/simple_plots_patch.py

#put the simulation into the calsim postgres database
python $pp2/sim2db.py

#make some plots of the work so far without color correcting
python $pp2/do_plots.py 

#save the plots
rm -r Iter0
mkdir Iter0
mv *.png *.eps delta* Iter0

#Make a color-correction map, apply the map to generate star_obs_corrected.dat, ses the current directory name as the name of the database
python $pp2/do_colorc.py

#rerun the photometric calibration using the corrected mags
${SELFCAL_SOLVER_DIR}/bin/simple.x star_obs_corrected.dat test 1.e-7 10000

#put the new test_* files into the db and make some plots
python $pp2/simple_plots_star.py
python $pp2/plots_fromdb.py

#save the first round plots
rm -r Iter1
mkdir Iter1
mv *.png *.eps delta* Iter1

#Do another round of color map generation, application and plot making
python $pp2/do_colorc.py

#rerun the photometric calibration using the corrected mags
${SELFCAL_SOLVER_DIR}/bin/simple.x star_obs_corrected.dat test 1.e-7 10000

python $pp2/simple_plots_star.py
python $pp2/plots_fromdb.py

##lather, rinse repeat until you think the fit converges.
rm -r Iter2
mkdir Iter2
mv *.png *.eps delta* Iter2

#Do another round of color map generation, application and plot making
python $pp2/do_colorc.py

#rerun the photometric calibration using the corrected mags
${SELFCAL_SOLVER_DIR}/bin/simple.x star_obs_corrected.dat test 1.e-7 10000

python $pp2/simple_plots_star.py
python $pp2/plots_fromdb.py
