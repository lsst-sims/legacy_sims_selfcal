#               Module : MJDtools
# Conversions    Modified Julian Day <-> Gregorian [Year, Month, Day, dayfract]
#                                     ->frommjd
#                                     <-tomjd
# includes Day fraction to fifth decimal
#  based on fromordinal and toordinal in datetime
from datetime import *
from numpy import mod
#_____________________________________________________________
#                frommjd
def frommjd(mjd):  
        gdat = date.fromordinal(int(mjd)+678576) # add ordinal of MJD origin
        gymdf = [int(str(gdat).split('-')[0]),\
                 int(str(gdat).split('-')[1]),\
                 int(str(gdat).split('-')[2]),\
                 mod(mjd,1)]
        return gymdf # [ year, month, day, dayfraction]
#_______________________________________________________________
#                tomjd
def tomjd(gymdf):  # input must be [year,month,day,dayfract]
        gdat = date(gymdf[0],gymdf[1],gymdf[2]).toordinal()
        mjdat = gdat - 678576
        mjd = mjdat+gymdf[3]
        return mjd
#______________________________________________________________
#                 ssn
# Given argument mjd, ssn returns in season
#                                - [0]current season number (0,1,2,3)
#                                - [1]season age at mjd (days)
#                                - [2] season length
#                                
def ssn(mjd):
        cur_year = frommjd(mjd)[0]
        if (mjd - tomjd([cur_year,12,21,0.0])) > 0 :
                cur_year = cur_year + 1
        s0 = tomjd([cur_year-1,12,21,0.0])
        s1 = tomjd([cur_year,3,21,0.0]) 
        s2 = tomjd([cur_year,6,21,0.0])
        s3 = tomjd([cur_year,9,21,0.0])
        s4 = tomjd([cur_year,12,21,0.0])
        sld  = [s1-s0, s2-s1, s3-s2, s4-s3] # season duration
        # derive  ssn = [ season at mjd, season age] at mjd
        if (mjd-s4)>0:
                ssn=[0,mjd-s4-1]   #   Winter (North hem.) Summer (South)
        elif (mjd-s3)>0 :
                ssn=[3,mjd-s3-1]   #   Spring (N) Fall (S)
        elif (mjd-s2)>0 :
                ssn=[2,mjd-s2-1]   #   Summer (N) Winter (S)
        elif (mjd-s1)>0 :
                ssn=[1,mjd-s1-1]   #   Fall (N)   Spring (S)
        else :
                ssn=[0,mjd-s0-1]   #   Winter (N) Summer (S)
        # add season length as apropriate --------------
        season= ssn + [sld[ssn[0]]]
        return season
#------------------------------------------------------
#                    sns
# Derive sequence of seasons, to occure between mjdstart and mjdend
# given season age and season length at MJDS
# 
#
# call as :  seasons = ssns(MJDstart,MJDend)

def ssns(MJDS,MJDE):
        season = ssn(MJDS)
        seas = [season]
        mjdsend = MJDS + season[2] - season[1]
        while mjdsend < MJDE :
                season = ssn(mjdsend)
                seas = seas + [season]
                mjdsend = mjdsend + season[2]
        return seas

