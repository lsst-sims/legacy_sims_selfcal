#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#"""                                                                   """
#"""                shape_store_spectra.py                                                   """
#"""               based on dmag_sequence_V2.py :                      """
#"""                                                                   """
#       compact and reformat .plt files
#     introduce reference spectrum

#--------------------------------------------------------------------------
# Read :
#   
#       - reference transmission spectrum at zenith
#
#       - Modtran atmosphere spectra for the sequence of pointings
#   
# 
#
# Write spectra columnwize in .csv## blockfiles of 1001 spectra
#  first column is wavelength, second is some reference spectrum 
# a summary of blockfile properties  in the csp_descript file
# 
#=========================================================================

from numpy import *

localpath =  '/groups/LSST/creze/Calsimwork/'
mlsstpath = localpath +'mkatmos/SpData/'

expname = 'Opsim3.61.02sel'

specpath = mlsstpath + expname + '/'
runpath = mlsstpath + 'Runparm/'
#======================================================================
#
#              TRANSMISSION SPECTRUM DATA ACQUISITION
#
#----------------------------------------------------------------------

nsbin = 5212
#----------------
# Read  reference transmission spectrum at zenith (return as 'rsp')
#----------------
#
firef = mlsstpath+'Ref400/Ref400.plt'
refsp = open(firef,'r')
rsp = []
for lsp in refsp:
    if lsp[7]=='$': 
        refsp.seek(0,2)  # stop reading at end of first spectrum
    else :
        rsp = rsp + [lsp[15:]]# rsp is list of 5212 reference transmittances
refsp.close()
		
#---------
# Read spectrum data for sequence (from .plt file)
#---------
#
fispec = specpath + expname +'.plt'

#
trsp = open(fispec,'r')
transpec = list(trsp)
trsp.close()
nlines = len(transpec)
#
# Put bin wavelengths in Lamspec[:,0]
#     reference transmission case in Lamspec[:,1]
#     transmission case 1 in Lamspec[:,2]
#     transmission case 2 in Lamspec[:,3]
#     ....
i1=-1   # i1 will be wavelength index
i2 = 1  #  i2 will be spectrum index (1 is reference spectrum)
lignes = [""]*nsbin
headline = [""]*nsbin
#blocksize = 1000  ############
blocksize = 1000    
#
# nlines = 100*(nsbin+1) # limitation for test
#
nblock = 0
for i in range(0,nlines):        #
    lf = transpec[i]
    if lf[7]<>'$':  # no $ in column 8 means new wavelength next
        i1=i1+1
        if i2== 1:
            # store wavelengthes, reference spectrum and first spectrum
            headline[i1] = headline[i1] + lf[:15].strip()+'$ '\
                        + rsp[i1][:-1].strip()+' '#+'$ '+ lf[15:].strip()
        if i2 > 0:# add transmission of new spectrum at wv + $ separator
            if mod(i2,blocksize) == 1:
                lignes[i1]=headline[i1]
            lignes[i1] =lignes[i1][:-1]+'$ '+ lf[15:].strip()
    if lf[7]=='$':  #  $ in column 8 means new spectrum next
        if mod(i2,blocksize)==0 or i== nlines-1 :
           nblock=nblock+1
           filename = specpath + expname + 'csp' + str(nblock)+'.dat'
           filout=open(filename,'w')
           for line in lignes :
               line = line + '\n'
               filout.write(line)           
           filout.close()
        i2 = i2+ 1
        # the last i2 gives the number of spectra (inlcuding reference)
        i1 = -1   
#
nspec = i2
lastblock = len(lignes[1].split('$'))
filename = specpath + expname + 'csp_control.dat'
contlist = ['sequence name = '+ str(expname)+'\n']
contlist = contlist + ['blocksize = '+ str(blocksize)+'\n']
contlist = contlist + ['totalspectra = '+ str(i2)+'\n']
contlist = contlist + ['bins = '+ str(nsbin)+'\n']
contlist = contlist + ['block number = '+ str(nblock)+'\n']
contlist = contlist + ['spectra per block = '+ str(blocksize)+'\n']
contlist = contlist + ['last block  = '+ str(lastblock-2)+'\n']
filout=open(filename,'w')
for txt in contlist :
    filout.write(txt)
filout.close()

#


