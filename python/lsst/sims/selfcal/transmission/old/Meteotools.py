#============================================================================
#
#               Library : Meteotools
#
# This tool box provides basic tools to simulate random   
#  meteorological series, currently (01-2009) only modules mkseg and mtsearch
#  are used by MakeAtmos to create and retrieve random time segments.
#  The remainder are attempts to implement continuity
#  of the meteorological variations up to the first derivative while
#  only zero order continuity is obtained.
#

from random import *
from numpy import *

#============================================================================
#
#                        mkseg
#
#       Make list of random exponential pseudo periods (segments)
#       inside interval [tstart tend]
#     for oscillating phenomena, there is the possibility that
#     semi periods up and semi periods down may be given different typical
#     durations ( ttypup /ttypdown )
#
#    Input  argument is tlims = [[tstart, tend], [ttypdown,ttypup]]
#       
#    Output is the list of intervals, each caracterized by ;
#       [t_cur, t_len,iupdown,afact,lastafact]
#       - t_cur is the start time of the current segment
#       - t_len is the time length of the current segment
#               t_len is random exponential with caracteristic scale
#               ttypup if idx = 1 ; ttypdown if idx = 0
#       - iupdown says whether variation is up or  down during segment
#       - afact is a (non negative, dimensionless)
#               random modulation of the parameter amplitude
#       - lastafact is the integrated memory of previous excursions
#               in amplitude, it is used to secure continuity across      
#               segment transition. This complexity is due
#               to the fact that segment transitions are operated
#               at extrema
#
#===========================================================================
#
def mkseg(tlims):
        tstart = tlims[0][0]
        tend = tlims[0][1]
        tcur = tstart      # first segment start
        lastafact = tlims[2][0]       
        iupdown = tlims[2][1]
        mtlis = []
        #print 'called mkseg between tstart ',tstart,' and tend ',tend 
        while tcur < tend :
                iupdown = -iupdown # switch current iupdown
                #print 'iupdown = ',iupdown
                print 'tcur = ',tcur
                idx = int( 1 + iupdown / 2. ) # idx is 0 (down) or 1 (up)
                tsca = tlims[1][idx]
                tlen = expovariate(1./tsca) # seg. length random expon
                tlen = max([tlen,0.3*tsca])
                tlen = min([tlen,3.*tsca])
                afact = max(0.1, normalvariate(1.,0.5)) # normal random amplitude factor
                #  
                #afact = 0.5 # dummy value for tests
                #
                mtlis = mtlis + [[tcur,tlen,iupdown,afact,lastafact]] #add one segment to list
                tcur = tcur + tlen # next segment start is end of current
                lastafact = lastafact - iupdown * afact # integrated memory of
                                   # previous amplitude excursions needed
                                   # to compute behaviour in next segment 
        return mtlis            # return list of segments
#________________________________________________________________
#
#                      mtsearch
#
#     Search in list of intervals, segment where tt is located
#     Input argument is ttmtl [tt,[mtlist]]where mtlist is list of segments
#      as created 
#      output  [tstart,tlength,idx]
#
def mtsearch(ttmtl):
        tt = ttmtl[0]
        mtlist  = ttmtl[1]
        ii = 0
        while mtlist[ii][0] < tt >= mtlist[ii][0]+mtlist[ii][1]:
                ii = ii + 1
                if ii == len(mtlist):
                        tt = mtlist[ii][0]+mtlist[ii][1]
        segment = mtlist[ii]
        return segment
#_______________________________________________________MJDE_________
#
#                      mtsmooth
#
#     Same as mtsearch 
#     But returns in addition previous and following segments
#     Output : 
#    [tstart0,tlength0,idx0],[tstart1,tlength1,idx1],[tstart2,tlength2,idx2]]
#       where index 1 refers to current segment
def mtssmooth(ttmtl):
        tt = ttmtl[0]
        mtlist  = ttmtl[1]
        ii = 0
        while mtlist[ii][0] < tt >= mtlist[ii][0]+mtlist[ii][1]:
                ii = ii + 1
        i0 = ii-1
        if i0<0 :
                seg0 = [mtlist[0][0]-mtlist[0][1]] + mtlist[0][1:]
        else :
                seg0 = mtlist[i0]
        i2 = ii+1
        if i2>=len(mtlist) :
                seg2 = [mtlist[ii][0]-mtlist[ii][1]] + mtlist[ii][1:]
        else :
                seg2 = mtlist[i2]
        segments = [seg0,mtlist[ii],seg2]
        return segments

#________________________________________________________________
#
#                      sinexp 
#
#       Computes  sine variations at regular absissae with period
#       changing over a series of intervals inside tstart tend
#       Input is arg = [ nabs , tlims ]
#        where tlims = tlims = [[tstart,tend],[[ttypup,ttypdown]]
#        typical test value : tlims = [[0.,20.],[2.,5.]]
#        calls mkexp ones to generate list of random exponential segments
#

def sinexp(arg):
        nabs = arg[0]
        tt=range(nabs)
        tlims = arg[1]
        tstep = (tlims[0][1]-tlims[0][0]) / nabs
        tt=array(tt) * tstep
        yy = zeros(len(tt))
        mtl = mkexp(tlims)
        #
        # for jj in range(0,len(mtl)):
        #       print jj , mtl[jj]
        #
        for ii in range(0,len(tt)):
                ttmtl = [tt[ii],mtl]
                seg = mtsearch(ttmtl)
                phase = (seg[2]+(tt[ii]-seg[0])/seg[1])*pi
                yy[ii] = sin(phase)
        curv=[tt,yy]
        return curv

#________________________________________________________________
#
#                      sinexps 
#
#      same as sinexp but phase is smoothed so as to get
#       continuity of first derivative through segment transition
#       Input is arg = [ nabs , tlims ]
#        where tlims = tlims = [[tstart,tend],[[ttypup,ttypdown]]
#        typical test value : tlims = [[0.,20.],[2.,5.]]
#        calls mkexp ones to generate list of random exponential segments
#

def sinexps(arg):
        nabs = arg[0]
        tt=range(nabs)
        tlims = arg[1]
        tstep = (tlims[0][1]-tlims[0][0]) / nabs
        tt=array(tt) * tstep
        yy = zeros(len(tt))
        mtl = mkexp(tlims)
        #
        # for jj in range(0,len(mtl)):
        #       print jj , mtl[jj]
        #
        for ii in range(0,len(tt)):
                ttmtl = [tt[ii],mtl]
                segs = mtsmooth(ttmtl)
                
                phase = (seg[2]+(tt[ii]-seg[0])/seg[1])*pi
                yy[ii] = sin(phase)
        curv=[tt,yy]
        return curv
