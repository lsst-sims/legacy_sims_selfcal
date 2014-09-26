# -*- coding: cp1252 -*-
#-------------------------------------
# module  Statphot.py
#---------------------------------------
# Version established  on March. 21 2010
#         finalized   on March 25 2010
# Purpose :
#       Visualize and Stat phtometric results of a sequence of MODTRAN
#        simulations.
# 
# Input :
#         - sequence of parameters used as input to Modtran
#                   from file  Runpar/..._parmlist.dat
#         - photometric data from a sequence of modtran simulations
#                  from file   ...._dmag.dat
#                      as produced by program dmag_SequenceV5.dat
# 
# Output:
#       curves of dmag's versus time compared to parms versus time
#       plots of dmag's versus / parms
#       histograms, mean, std, min, max of dmags 
#           for various parm selections
#       There is a handle under the # parameter selection criterion
#       which allows to activate a selection of low airmass cases
#------------------------------------------------------------
#
import numpy
import pylab
import Astrotools
from mysort import csort

# common items
localpath =  '/groups/LSST/creze/Calsimwork/'
mlsstpath = localpath +'mkatmos/SpData/'




#----------
# experience set name and description file identification
#      THIS parameter MUST BE CHANGED for each new experience
#
seqname = 'Opsim3.61.02sel'
seqpath = mlsstpath + seqname +'/'
runpath = mlsstpath + 'Runparm/'


sortbyairmass = False
d2r=numpy.pi/180.

#
#==================================================================
#                 functions
#

#------------------------------------------------------------------
#
#            histogram plot routine
#

def subhplot(hplot):  # create histogram plots
    hh = hplot[0]
    xmin = hplot[1][0]
    xmax = hplot[1][1]
    xtipos = hplot[1][2]
    ymin = hplot[2][0]
    ymaxf = hplot[2][1]
    fmx = hplot[2][2]
    ytiposf = hplot[2][3]
    hname = hplot[3]    
    sbp = hplot[4]
    nsub = numpy.mod(sbp,10)
    pylab.plot(hh[1],hh[0],linestyle='None',linewidth=0.1 )
    yyl = numpy.log10(max(hh[0]))
    ylf = numpy.fix(yyl)
    ymax = (numpy.round(10**(yyl - ylf)))*(10**ylf)
    ystep=numpy.round(ymax-ymin)/4
    ytipos= numpy.arange(ystep,4*ystep,ystep)
    if fmx == True :
        ymax = ymaxf
        ytipos = ytiposf
    pylab.text(xmin+0.1*(xmax-xmin),0.7*ymax,hname)
    pylab.axis([xmin,xmax,ymin,ymax],linewidth=0.1)
    #if numpy.mod(nsub,2)==1 :
    if numpy.mod(nsub,2)<2 :
        pylab.yticks(ytipos,weight='light',fontsize=9)
    if nsub >= 5:
        pylab.xticks(xtipos,weight='light',fontsize=9)
    xtiblank=(((''),)*len(xtipos))
    ytiblank=(((''),)*len(ytipos))
    #if numpy.mod(nsub,2)==0 :
    #    pylab.yticks(ytipos,ytiblank)
    if nsub < 5:
        pylab.xticks(xtipos,xtiblank)
    return

#------------------------------------------------------------------------
#
#                   x,y plot routine
#

def subxyplot(xyplot):  # create x,y plots
    xy = xyplot[0]
    xmin = xyplot[1]
    xmax = xyplot[2]
    ymin = xyplot[3]
    ymax = xyplot[4]
    ytipos = xyplot[5]
    pname = xyplot[6]
    sbp = xyplot[7]
    nsub = numpy.mod(sbp,10)
    pylab.plot(xy[0],xy[1],'.b',linestyle='None',markersize=1)
    pylab.axis([xmin,xmax,ymin,ymax],linewidth=0.1)
    xstep=(xmax-xmin)/5
    ystep= (ymax-ymin)/5
    xtipos= numpy.arange((xmin+xstep),(xmin+5*xstep),xstep)
    #ytipos= numpy.arange(0.,3*ystep,ystep)
    pylab.text(xmin+(xstep/2.),0.7*ymax,pname)
    #if numpy.mod(nsub,2)==1 :#to be set if vscale not to be repeted in column 2
    if numpy.mod(nsub,2)<2 :
        pylab.yticks(ytipos,weight='light',fontsize=9)
    if nsub >= 5:
        pylab.xticks(xtipos,weight='light',fontsize=9)
    xtiblank=(((''),)*len(xtipos))
    ytiblank=(((''),)*len(ytipos))
    #if numpy.mod(nsub,2)==0 :
    #    pylab.yticks(ytipos,ytiblank)
    if nsub < 5:
        pylab.xticks(xtipos,xtiblank)
    return



#===============================================================
#
#              icol
#
#  find column rank in dmag file for photometry pass band
#         identified by filter name
#  
#  
#
# 
def icol(name) :
    dmag_columnhead = ['u','g','r','i','z','y']  
    lnam = dmag_columnhead
    icolumn=-1
    for cn in lnam :
        icolumn = icolumn +1
        #print icolumn, '  name',name,'  colname',cn,
        if name == cn :
            return icolumn

#=======================================================================
#
#           colquery
#
# colquery is same as icol, but list of columnheads is passed as query[1]
def colquery(query):
    name = query[0]
    lnam = query[1]
    icolumn=-1
    for cn in lnam :
        icolumn = icolumn +1
        if name == cn :
            return icolumn
#
#
#=========================================================
#        Input acquisition
#----------------------------------------------------------------

#-------
# parameter  selection criterio
Nosel = False
Nosel = True
if Nosel == True :
    selsuffix = ''
if Nosel ==False :
    parsel = 'ZANGLE'
    selval=30.
    seltest = 'le'# mind : Test must be set accordingly in Parameter sel section 
    selsuffix = parsel[:2]+'.'+seltest+'.'+str(int(selval))
    selval = float(selsuffix[-2:])

#
#----------          
#  photometric data acquisition (returned in list lphp and array phot)
photfile = seqpath + seqname+'_dmag.dat'
fphp = open(photfile,'r')
lphp = list(fphp)[4:]  # switch file content to list and skip first 5 lines
fphp.close()
nruns = len(lphp)
#   extract photometry columns data in array "phot"
ndat = len(lphp[0][14:-2].split('|'))              
irun = -1
phot = [[0.]*ndat]*nruns
ident = [('')]*nruns
for line in lphp[1:] :
    fldat = []
    if line[:12].strip()<> 'Reference': # exlude reference lines from data
        irun = irun + 1
        ldat = line[14:-2].split()
        for ii in range(0,len(ldat)):
            fldat = fldat + [float(ldat[ii])]
        phot[irun]= fldat
        ident[irun]= line.split(':')[0].strip()
nrun = irun + 1
phot=numpy.array(phot[:nrun])
ident = ident[:nrun]

#----------          
#  parameter data acquistion (returned in line list lparl and array parma)
parmfile = runpath + seqname + '_parmlist.dat'
fpar = open(parmfile,'r')
lparl = list(fpar)
fpar.close()
nparun = len(lparl)
if nparun <> len(ident):
    print '***WARNING*** lengthes of parmlist and photlist missmatch'

# names of parms to be extracted
parnaml = ['expMJD', 'ZANGLE','MODEL','O3','H2O','ISEAS','ICLD','VIS']
nparm = len(parnaml)
parma = [[0.]*nparm]*nrun
parma = numpy.array(parma)
irun=-1

for line in lparl :
    irun = irun + 1
    parid = line.split('$')[0].strip()
    if parid <> ident[irun]:
        print '***WARNING*** in trouble with identifiers'\
              , irun ,' parid ',parid,' ident ',ident
    lrpar = line.split('$')
    for ipa in range(0,nparm):
        parname = parnaml[ipa]
        for parmstr  in lrpar[1:] :
            pa = parmstr.split('=')
            if pa[0].strip()==parname :
               parma[irun,ipa]=float(pa[1].strip())    
#
#                      ....
#
if sortbyairmass == True :
    phot = numpy.transpose(csort((numpy.transpose(phot),icol('Zangle'))))

# photometric indices
u_g=phot[:,icol('u')]-phot[:,icol('g')]
g_r=phot[:,icol('g')]-phot[:,icol('r')]
g_i=phot[:,icol('g')]-phot[:,icol('i')]
r_i=phot[:,icol('r')]-phot[:,icol('i')]
r_z=phot[:,icol('r')]-phot[:,icol('z')]
r_y=phot[:,icol('r')]-phot[:,icol('y')]
colar= numpy.transpose([u_g,g_r,g_i,r_i,r_z,r_y])
colnaml = ['u_g','g_r','g_i','r_i','r_z','r_y']
#===========================================================================
#
#           Parameter selection area
#

if Nosel == False :
    ipsel = colquery([parsel,parnaml])
    Vsel = parma[:,ipsel]
    r_colar = []
    r_phot = []
    r_parma = []
    isel=-1
    for idx in range(0,len(Vsel)):
        if Vsel[idx] <= selval :
            r_colar = r_colar + [colar[idx]]
            r_phot = r_phot + [phot[idx]]
            r_parma = r_parma + [parma[idx]]
    colar =numpy.array( r_colar)
    phot = numpy.array(r_phot)
    parma = numpy.array(r_parma)

#===========================================================================
#
#                Plot area
#
#=================================
#               1. Histograms
#
#----------------
# 1.1     global, dmags
pylab.figure(1)
pylab.clf()
bins= numpy.arange(-0.1,0.35,0.02)
xmax =0.3
xmin = -0.1
xtipos = [0.0,0.1,0.2]
#
ymin = 0
ymaxf = 1500
ytiposf = [250,500,750,1000,1250]
fmx = True # if True ymax will be forced to ymf  
#
lisnam= ['u','g','r','i','z','y']
sbp0= 321
for i in range(0,len(lisnam)) :
    hname = lisnam[i]
    sbp=sbp0+i
    pylab.subplot(sbp)
    hh = pylab.hist(phot[:,icol(hname)],bins,linewidth=0.1)
    hplot = (hh,[xmin,xmax,xtipos],[ymin,ymaxf,fmx,ytiposf],hname,sbp)
    subhplot(hplot)
yyl = numpy.log10(max(hh[0]))
ylf = numpy.fix(yyl)
ymax = (numpy.round(10**(yyl - ylf)))*(10**ylf)
if fmx == True :
    ymax = ymaxf
pylab.text((xmin-1.5*xmax),3.5*ymax,'Histogram of atmosphere generated magnitude '\
           +'deviations (global) : '+seqname+'  '+selsuffix, fontsize = 9)
hfig = seqpath + seqname+selsuffix + '_dmag_hist.eps'
pylab.savefig(hfig)
#
#
#----------------
# 1.2   global (colours)
pylab.clf()
bins= numpy.arange(-0.1,0.15,0.01)
xmax =0.15
xmin = -0.10
xtipos = [-0.05,0.0,0.05,0.10]
#
ymin = 0
ymaxf = 1500
ytiposf = [250,500,750,1000,1250]
fmx = True # if True ymax will be forced to ymf
#
sbp0= 321
i=-1
for i in range(0,len(colnaml)) :
    hname = [colnaml[i],colnaml]
    sbp=sbp0+i
    pylab.subplot(sbp)
    hh = pylab.hist(colar[:,colquery(hname)],bins,linewidth=0.1)
    hplot = (hh,[xmin,xmax,xtipos],[ymin,ymaxf,fmx,ytiposf],hname[0],sbp)
    subhplot(hplot)
yyl = numpy.log10(max(hh[0]))
ylf = numpy.fix(yyl)
ymax = (numpy.round(10**(yyl - ylf)))*(10**ylf)
if fmx == True :
    ymax = ymaxf
pylab.text((xmin-2.*xmax),3.5*ymax,'Histogram of atmosphere generated colour '\
           +'deviations (global) : '+seqname+'  '+selsuffix, fontsize = 9)
hfig = seqpath + seqname+selsuffix + '_dcol_hist.eps'
pylab.savefig(hfig)

#----------------
# 1.3   selected by IVULC
#
#----------------
# 1.4   selected by Zangle
#

#=================================
#               2. magnitude plots
#
ymin= - 0.1
ymax =  0.4
ytipos = [0.0,0.1,0.2,0.3]
ytit = 1.65
#----------------
# 2.1     versus time
#
# pylab.figure(2)
pylab.clf()
xparn = 'expMJD'
ipa = colquery([xparn,parnaml])
xval = parma[:,ipa]
xval = xval - numpy.fix(numpy.min(xval))
xmin = 0
xmax = numpy.fix(numpy.max(xval)+10. - numpy.mod(numpy.max(xval),10))
sbp0 = 321
for i in range(0,len(lisnam)) :
    hname = lisnam[i]
    sbp=sbp0+i
    pylab.subplot(sbp)
    yval = phot[:,icol(hname)]
    xy = [xval,yval]
    xyplot = (xy,xmin,xmax,ymin,ymax,ytipos,hname,sbp)
    subxyplot(xyplot)
pylab.text(-1300,ytit,' Atmosphere generated magnitude '\
           +'deviates vs time (days since mjd 49353): '\
           + seqname+'  '+selsuffix , fontsize = 9)
hfig = seqpath + seqname+selsuffix + '_dmag_time.eps'
pylab.savefig(hfig)


#----------------
# 2.2     versus H2O
# 
#
pylab.clf()
xparn = 'H2O'
ipa = colquery([xparn,parnaml])
xval = parma[:,ipa]
xval = xval - numpy.fix(numpy.min(xval))
xmin = 0.4
xmax = 1.6
sbp0 = 321
for i in range(0,len(lisnam)) :
    hname = lisnam[i]
    sbp=sbp0+i
    pylab.subplot(sbp)
    yval = phot[:,icol(hname)]
    xy = [xval,yval]
    xyplot = (xy,xmin,xmax,ymin,ymax,ytipos,hname,sbp)
    subxyplot(xyplot)
pylab.text(-1.,ytit,' Atmosphere generated magnitude '\
           +'deviates vs H2O : '\
           + seqname+'  '+selsuffix , fontsize = 9)
hfig = seqpath + seqname+selsuffix + '_dmag_H2O.eps'
pylab.savefig(hfig)


#----------------
# 2.3     versus VIS
#
# 
# 
# 
pylab.clf()
xparn = 'VIS'
ipa = colquery([xparn,parnaml])
xval = parma[:,ipa]
xval = xval - numpy.fix(numpy.min(xval))
xmin = 0.
xmax = 60.
sbp0 = 321
for i in range(0,len(lisnam)) :
    hname = lisnam[i]
    sbp=sbp0+i
    pylab.subplot(sbp)
    yval = phot[:,icol(hname)]
    xy = [xval,yval]
    xyplot = (xy,xmin,xmax,ymin,ymax,ytipos,hname,sbp)
    subxyplot(xyplot)
pylab.text(-50.,ytit,' Atmosphere generated magnitude '\
           +'deviates vs VIS : '\
           + seqname+'  '+selsuffix , fontsize = 9)
hfig = seqpath + seqname+selsuffix + '_dmag_VIS.eps'
pylab.savefig(hfig)
#----------------
# 2.4     versus Zangle
#
# 
# 
# 
pylab.clf()
xparn = 'ZANGLE'
ipa = colquery([xparn,parnaml])
xval = parma[:,ipa]
xval = xval - numpy.fix(numpy.min(xval))
xmin = 0.
xmax = 75.
sbp0 = 321
for i in range(0,len(lisnam)) :
    hname = lisnam[i]
    sbp=sbp0+i
    pylab.subplot(sbp)
    yval = phot[:,icol(hname)]
    xy = [xval,yval]
    xyplot = (xy,xmin,xmax,ymin,ymax,ytipos,hname,sbp)
    subxyplot(xyplot)
pylab.text(-xmax,ytit,' Atmosphere generated magnitude '\
           +'deviates vs Zangle : '\
           + seqname+'  '+selsuffix , fontsize = 9)
hfig = seqpath + seqname+selsuffix + '_dmag_Zangle.eps'
pylab.savefig(hfig)

#=================================
#               3. colour plots
#

ymin= - 0.1
ymax =  0.2
ytipos = [-0.05,0.0,0.05,0.10,0.15]
ytit = 0.95
#----------------
# 3.1     versus time
#
# pylab.figure(2)
pylab.clf()
xparn = 'expMJD'
ipa = colquery([xparn,parnaml])
xval = parma[:,ipa]
xval = xval - numpy.fix(numpy.min(xval))
xmin = 0
xmax = numpy.fix(numpy.max(xval)+10. - numpy.mod(numpy.max(xval),10))
sbp0 = 321
for i in range(0,len(colnaml)) :
    hname = [colnaml[i],colnaml]
    sbp=sbp0+i
    pylab.subplot(sbp)
    yval = colar[:,colquery(hname)]
    xy = [xval,yval]
    xyplot = (xy,xmin,xmax,ymin,ymax,ytipos,hname[0],sbp)
    subxyplot(xyplot)
pylab.text(2.*xmin-1.1*xmax,ytit,' Atmosphere generated colour '\
           +'deviates vs time (days since mjd 49353): '\
           + seqname+'  '+selsuffix , fontsize = 9)
hfig = seqpath + seqname+selsuffix + '_dcol_time.eps'
pylab.savefig(hfig)


#----------------
# 3.2     versus H2O
# 
#
pylab.clf()
xparn = 'H2O'
ipa = colquery([xparn,parnaml])
xval = parma[:,ipa]
xval = xval - numpy.fix(numpy.min(xval))
xmin = 0.4
xmax = 1.6
ymin= - 0.1
ymax =  0.2
sbp0 = 321
for i in range(0,len(colnaml)) :
    hname = [colnaml[i],colnaml]
    sbp=sbp0+i
    pylab.subplot(sbp)
    yval = colar[:,colquery(hname)]
    xy = [xval,yval]
    xyplot = (xy,xmin,xmax,ymin,ymax,ytipos,hname[0],sbp)
    subxyplot(xyplot)
pylab.text(2.*xmin-1.1*xmax,ytit,' Atmosphere generated colour '\
           +'deviates vs H2O : '\
           + seqname+'  '+selsuffix , fontsize = 9)
hfig = seqpath + seqname+selsuffix + '_dcol_H2O.eps'
pylab.savefig(hfig)


#----------------
# 3.3     versus VIS
#
# 
# 
# 
pylab.clf()
xparn = 'VIS'
ipa = colquery([xparn,parnaml])
xval = parma[:,ipa]
xval = xval - numpy.fix(numpy.min(xval))
xmin = 0.
xmax = 60.
ymin= - 0.1
ymax =  0.2
sbp0 = 321
for i in range(0,len(colnaml)) :
    hname = [colnaml[i],colnaml]
    sbp=sbp0+i
    pylab.subplot(sbp)
    yval = colar[:,colquery(hname)]
    xy = [xval,yval]
    xyplot = (xy,xmin,xmax,ymin,ymax,ytipos,hname[0],sbp)
    subxyplot(xyplot)
pylab.text(2.*xmin-1.1*xmax,ytit,' Atmosphere generated colour '\
           +'deviates vs VIS : '\
           + seqname+'  '+selsuffix , fontsize = 9)
hfig = seqpath + seqname+selsuffix + '_dcol_VIS.eps'
pylab.savefig(hfig)
#----------------
# 3.4     versus Zangle
#
# 
# 
# 
pylab.clf()
xparn = 'ZANGLE'
ipa = colquery([xparn,parnaml])
xval = parma[:,ipa]
xval = xval - numpy.fix(numpy.min(xval))
xmin = 0.
xmax = 75.
ymin= - 0.1
ymax =  0.2
sbp0 = 321
for i in range(0,len(colnaml)) :
    hname = [colnaml[i],colnaml]
    sbp=sbp0+i
    pylab.subplot(sbp)
    yval = colar[:,colquery(hname)]
    xy = [xval,yval]
    xyplot = (xy,xmin,xmax,ymin,ymax,ytipos,hname[0],sbp)
    subxyplot(xyplot)
pylab.text(2.*xmin-1.1*xmax,ytit,' Atmosphere generated colour '\
           +'deviates vs Zangle : '\
           + seqname+'  '+selsuffix , fontsize = 9)
hfig = seqpath + seqname+selsuffix + '_dcol_Zangle.eps'
pylab.savefig(hfig)
#----------------
# 3.5     versus O3
#
# 
# 
# 
pylab.clf()
xparn = 'O3'
ipa = colquery([xparn,parnaml])
xval = parma[:,ipa]
xval = xval - numpy.fix(numpy.min(xval))
xmin = 0.5
xmax = 1.5
ymin= - 0.1
ymax =  0.2
sbp0 = 321
for i in range(0,len(colnaml)) :
    hname = [colnaml[i],colnaml]
    sbp=sbp0+i
    pylab.subplot(sbp)
    yval = colar[:,colquery(hname)]
    xy = [xval,yval]
    xyplot = (xy,xmin,xmax,ymin,ymax,ytipos,hname[0],sbp)
    subxyplot(xyplot)
pylab.text(2.*xmin-1.1*xmax,ytit,' Atmosphere generated colour '\
           +'deviates vs O3 : '\
           + seqname+selsuffix , fontsize = 9)
hfig = seqpath + seqname+selsuffix + '_dcol_O3.eps'
pylab.savefig(hfig)


