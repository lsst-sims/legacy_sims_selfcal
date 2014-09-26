# -*- coding: cp1252 -*-
#-------------------------------------
# module  Makecard.py
#---------------------------------------
# Version established and tested on Jan. 15 2009
#
# Read observing target inputs from xxx_parmlist.dat
# this .dat file must specify target airmass and time
# and every parameter requested by MODTRAN to fix the atmospheric conditions
# at the time of observations
# This Version includes the possibility to insert (or not) optional Card 2A
# Requested when cirrus clouds are present
# It also includes Card 2B to specify aerosol features
# The template in use must include the 2 cards 2A+ which place
#   the aerosol layers in altitude, this feature is essential
#   to force ground generated aerosols at the altitude of the telescope
# 
# 
#------------------------------------------------------------
#
from numpy import *
from copy import deepcopy
import readcard
import readparms
#
#------------------------------------------------------------
#        Input acquisition
#------------------------------------------------------------
#  paths to data
#                    Give below absolute path to mkatmos in principle it is the only change
localpath =  '/Users/ljones/lsst_home/Darwin/sims/atmosphere/'
mlsstpath = localpath +'mkatmos/'
runpath = mlsstpath + 'SpData/Runparm/'
parpath = mlsstpath + 'Parmfiles/'

#
# Modtran input Card template acquisition
#
card0file = parpath +'Cardtemplate.dat'
cardtemp = readcard.rdc(card0file)
ncardt = len(cardtemp)# get number of cards in template
nc0 = ncardt - 1

#=========================================================
# experience set name and description file identification
#      THIS parameter MUST BE CHANGED
#-------------------
expname = 'Opsim3.61.02sel'

#         ==========
#
#
#
# observing target and atmosphere input parameter acquisition
#

parmfile = runpath +expname+ '_parmlist.dat'
print 'parmfile %s' %(parmfile)
runs = readcard.rdc(parmfile) # runs will gather parameter string for all runs

################################################
#              generate new card list
#  ncl = cardtemp # this is the zero change case (new = template)
#
nrun = len(runs)  # number of runs in current set


parms=[]
for irun in range(0,nrun):
    # get parameter to be changed
    pars = runs[irun]    # isolate parameter string defining one run
    psp = pars.split('$')# split in list
    parname = []         #initialize list of parameter names
    parval = []          # initialize list of parameter values
    npar = len(psp)-1    # number of parameters to be changed for this run
    for ipar in range(1,npar):
        par = psp[ipar].split('=')
        parname = parname+[par[0].strip()]
        parval = parval+[par[1].strip()]
        paripar = [parname, parval]
    
    parms = parms + [paripar]

# the first index in parms selects runs (from 0 to nrun-1)
# the second index switches between name [0] and value [1]
# the third index selects parameters (from 0 to npar-1)
# npar(irun) is obtained as len(parms[irun-1][1])
#__________


# get list of parameters authorized for change with suitable format
#  and location in modtran cards

print "here"
parcatf = parpath +'formparm.dat'
parmcat = readparms.rdp(parcatf)

# create modified list of cards
newcards=[]
for irun in range(0,nrun) :
    newc = deepcopy(cardtemp)  
    ir =irun
    npar = len(parms[ir][1])
    for ipar in range(0,npar):
        for icat in range(0,len(parmcat)):
            if parms[ir][0][ipar]==parmcat[icat][0]:
                pnam = parmcat[icat][0]
                #print pnam
                pcard = parmcat[icat][1] # rank of card to be modified
                pform = parmcat[icat][2] # parameter format in card
                pstart = parmcat[icat][3]-1 # parameter location in card
                pend = parmcat[icat][4]     # end "      "        "  "
                # convert parameter value string to appropriate type
                if pform[-1]== 's':
                    parval = parms[ir][1][ipar]
                if pform[-1]== 'd':
                    parval = int(parms[ir][1][ipar])
                if pform[-1]== 'f':
                    parval = float(parms[ir][1][ipar])
                # substitute  parameter value at               
                newc[pcard] = newc[pcard][0:pstart]+\
                              pform %parval+\
                              newc[pcard][pend:]
    # set continuation card to 1 : yes (full set)

    newc[nc0]=newc[nc0][:4]+'1'+newc[nc0][5:]
    # print irun, 'newc', newc
    # print 'cardtemp', cardtemp
    newcards=newcards + [newc]
 
# set last continuation card to 0
newcards[nrun-1][nc0]= newc[nc0][:4]+'0'+newc[nc0][5:]
# store new card sequence in list ncl

ncl=[]
# define optional card2A and/or Card2B to be inserted after card2
# in case clouds are expected 
card2A =  ['   0.000   0.000   0.000      '\
          +'                              '\
          +'                    \n']#default values for Cirrus 18 or 19
# incase aerosols are non standard
card2B = ['    -1.000     2.000     1.000 '\
          +'                              '\
          +'                   \n']# aerosol fog near surface 
#                                           # decreasing with height
for irun in range(0,nrun):
    # insert card2A and/or Card2B whenever required
    c2= newcards[irun][3]
    ihaz = int(c2[2:5])
    icld = int(c2[20:26])
    ivsa = int(c2[25:30])
    adcard = 6   #this is the rank of the first extra card to be inserted
                 # should be 6 if card 2A+ is placed in the template set
    if ihaz > 0 and icld > 0 :
        newcards[irun] = newcards[irun][:adcard]+card2A+newcards[irun][adcard:]
        adcard = adcard + 1 
    if ihaz > 0 and ivsa > 0 :
        newcards[irun] = newcards[irun][:adcard]+card2B+newcards[irun][adcard:]       
    ncl = ncl + newcards[irun]
 
######################################################
#  write new card list in tnn.nn.tp5 file
newcardfil = open(mlsstpath+'SpData/'+expname+'/' + expname + '.tp5','w')
for ic in range(0,len(ncl)):
    newcardfil.write(str(ncl[ic]))

newcardfil.close()


    
               
    

