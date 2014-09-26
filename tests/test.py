import os
import pylab
import lsst.sims.atmosphere.transmission.modtranCards as mc
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass

import time
def dtime(time_prev):
   return (time.time() - time_prev, time.time())



def compare_line_by_line(file1='tmp.tp5', file2='ref_tmp.tp5'):
    # compare tmp.tp5 with tmp.tp5_ref
    file1 = open(file1, 'r')
    file2 = open(file2, 'r')
    for line1, line2 in zip(file1, file2):
        print 'new', line1.strip()
        print 'ref', line2.strip()
    file1.close()
    file2.close()
    return


# Instantiate object
m = mc.ModtranCards()

# Read default card template
m.readCardTemplate()
#print m._cardTemplate

# Read default parameter formats
m.readParameterFormats()
#print m._paramFormats

# Read a test input parameter file (in Michel's format)
paramList = m.readParamValues_M('tmp_parmlist.dat')

# Look at a few atmospheres
tests = 2

for i in range(tests):
    t = time.time()
    print "Generating atmosphere for a test case:"
    print paramList[i]
    
    # Write modtran input cards
    m.writeModtranCards(paramList[i], 'tmp')

    # Run modtran on this file.
    m.runModtran('tmp')
    
    dt, t = dtime(t)
    print "Generating atmosphere took %f seconds." %(dt)

    # Read the atmosphere back in and plot to screen. 
    atm = Bandpass()
#    atm.readThroughput('tmp.psc', wavelen_min=800, wavelen_max=1300, wavelen_step=0.001)
    atm.readThroughput('tmp.psc', wavelen_min=300, wavelen_max=1100, wavelen_step=0.1)
    pylab.plot(atm.wavelen, atm.sb, 'b-')
    pylab.xlabel('Wavelength (nm)')
    pylab.ylabel('Atmospheric Transmission')
    
    m.cleanModtran('tmp')


pylab.show()




