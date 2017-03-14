from __future__ import print_function
####################################################################
#  Written by: Peter Yoachim - UW - v1.0 10/21/10
#  Questions or comments, email : yoachim@uw.edu
#
#  
#  Python program to take the output of simSelfCalib and dump it into a
#  Postgres database.  The code can then query the DB to construct a
#  color correction map.
#
#  The code assumes you have a local postgres database named "calsim"
#  running.  
#
#
#
#    software requirements 
#         - numpy (version 1.2 or higher)
#         - pylab --trying to switch to plain matplotlib
#         - pgdb
#
####################################################################


import sys
import os
import time 
import mmap

import numpy as n
import matplotlib
matplotlib.use('Agg')
#import matplotlib
import pylab as pyl
#import pgdb as pg
import psycopg2 as pg

import lsst.sims.selfcal.analysis.useful_input as ui
import lsst.sims.selfcal.analysis.catalog_plots_utils as cpu
import lsst.sims.selfcal.analysis.catalog_plots as cp

pyl.rcParams['font.size']=16



def plot_contour_irregulargrid(x, y, z, nxbins=360, nybins=180, cb_title=None):
    """Make a contour plot of irregularly gridded data."""
    # Grid the data.
    dx = (n.max(x) - n.min(x))/float(nxbins)
    dy = (n.max(y) - n.min(y))/float(nybins)
    xi = n.arange(n.min(x), n.max(x), dx)
    yi = n.arange(n.min(y), n.max(y), dy)
    xi, yi = n.meshgrid(xi, yi)
    zi = pyl.griddata(x, y, z, xi, yi)
    # Contour the gridded data.
    CS = pyl.contour(xi, yi, zi)
    CS = pyl.contourf(xi, yi, zi, 200, cmap=pyl.cm.jet)
    cb = pyl.colorbar()
    cb.set_label(cb_title)
    # Plot data points
    #pyl.scatter(x, y, 'b.')
    return 

def sim2db(master_table_name, index_tag, masterfile='master_cal.dat', starfile='stardata.dat',
           bestfitfile='test_bestfit_Star.dat', bestfitpatchfile='test_bestfit_Patch.dat', 
           star_obs=None):
    """Take the simulation output from master_table.dat and the bestfit files. 
    Put them into a postgres database.  Index the tables on useful quanities (like StarID).  """

    bestfit_table_name = master_table_name+'_bestfit'
    bestfit_patch_table_name = master_table_name+'_bestfitpatch'

    #reformat simple.x output    
    os.system("perl -i -p -e 's/^\s+//;' "+bestfitfile) ##remove leading spaces
    os.system("perl -i -p -e 'tr/ //s;' "+bestfitfile)  ##remove  any extra spaces
    os.system("perl -i -p -e 's/^\s+//;' "+bestfitpatchfile) ##remove leading spaces
    os.system("perl -i -p -e 'tr/ //s;' "+bestfitpatchfile)  ##remove  any extra spaces

    #open the database
    print('opening connection')
    #db = pg.connect(dsn = 'localhost:calsim')
    db = pg.connect(host='localhost', dbname='calsim')
    cursor = db.cursor ()
    cwd=os.getcwd()

    #copy in the stardata
    if starfile != None:
        cursor.execute("DROP TABLE IF EXISTS "+master_table_name+"_stars")
        cmd = "CREATE TABLE "+master_table_name+"_stars (StarID int, StarDBID int, RA double precision, Dec double precision, StarTrueMag double precision, StarColorTrue double precision);"
        # This should be replaced with a "create table like" (and store table format in 'master')

        cursor.execute(cmd)
        cmd="COPY "+master_table_name+"_stars FROM '"+cwd+"/stardata.dat' DELIMITER ' ' HEADER CSV;"
        cursor.execute(cmd)
        cursor.execute('DROP INDEX IF EXISTS stardataid'+index_tag+' ;')
        cmd='CREATE INDEX stardataid'+index_tag+' on '+master_table_name+'_stars (StarID)'
        print(cmd)
        cursor.execute(cmd)
        db.commit()

    if masterfile != None:
        print('creating master table')
        cursor.execute ("DROP TABLE IF EXISTS "+master_table_name)
        print("DROP TABLE IF EXISTS "+master_table_name)
        cmd = "CREATE TABLE " + master_table_name + " ( \
        VisitID   int, \
        PatchID   int,\
        subPatchID  int,\
        StarID  int,\
        StarDBID  int,\
        StarObsMag  double precision, \
        StarMagErr  double precision, \
        dMag_var	 double precision, \
        dMag_snr	 double precision, \
        dMag_zp	 double precision, \
        dMag_zp_dist double precision, \
        dMag_color	  double precision, \
        dMag_color_dist  double precision, \
        dMag_sin	      double precision, \
        dMag_flat_sin	      double precision, \
        dMag_cloud            double precision, \
        dMag_cloud_image     double precision, \
        dMag_rr           double precision, \
        StarTrueMag      double precision, \
        StarTrueCol      double precision, \
        StarRA	      double precision, \
        StarDec	      double precision, \
        RAcen_fov	      double precision, \
        DecCen_fov	      double precision, \
        StarX	      double precision, \
        StarY	      double precision,\
        night         double precision,\
        dm5           double precision,\
        dMag_kep      double precision,\
        dMag_gainvar      double precision, \
        dmag_exptvar      double precision\
        );"

        print(cmd)
        cursor.execute(cmd)
        db.commit()

        ##copy the masterfile into the db
        cmd = "COPY " + master_table_name + " from '"+cwd+'/' \
              + masterfile + "' WITH DELIMITER ' ' HEADER CSV;"
        print('Putting the big master file in the DB')
        print(cmd)
        cursor.execute (cmd)
        db.commit()

        cursor.execute('DROP INDEX IF EXISTS starx_index'+index_tag+' ;')
        cursor.execute('DROP INDEX IF EXISTS stary_index'+index_tag+' ;')
        cursor.execute('DROP INDEX IF EXISTS color_index'+index_tag+' ;')
        cursor.execute('DROP INDEX IF EXISTS StarObsMag_index'+index_tag+' ;')
        cursor.execute('DROP INDEX IF EXISTS starid_index'+index_tag+' ;')
        cursor.execute('DROP INDEX IF EXISTS night_index'+index_tag+' ;')

        cmd = "CREATE INDEX starx_index" +index_tag +" on " +master_table_name + " (StarX);"
        cursor.execute (cmd)
        cmd = "CREATE INDEX stary_index" +index_tag +" on " +master_table_name + " (StarY);"
        cursor.execute (cmd)
        cmd = "CREATE INDEX color_index" +index_tag +" on " +master_table_name + " (StarTrueCol);"
        cursor.execute (cmd)
        cmd = "CREATE INDEX StarObsMag_index" +index_tag +" on " +master_table_name + " (StarObsMag);"
        cursor.execute (cmd)
        cmd = "CREATE INDEX starid" +index_tag +" on " +master_table_name + " (StarID);"
        cursor.execute (cmd)
        cmd = "CREATE INDEX night_index" +index_tag +" on " +master_table_name + " (night);"
        cursor.execute (cmd)

        cursor.execute('DROP INDEX IF EXISTS patchid_patch'+index_tag+' ;')
        cmd = "CREATE INDEX master_patchid" + index_tag +" on " +  master_table_name + " (PatchID);"
        cursor.execute (cmd)

        #XXXXXXX---make a ditty here that copies over the star_obs.dat file and then remake the index.
    
    #Create the bestfit stars table
    cursor.execute ("DROP TABLE IF EXISTS "+bestfit_table_name)
    print("DROP TABLE IF EXISTS "+bestfit_table_name)
    cmd = "CREATE TABLE " + bestfit_table_name + " (StarID int, StarMagCalib double precision);"
    print(cmd)
    cursor.execute (cmd)
    db.commit()

    ##Copy the bestfit stars stuff in
    cmd = "COPY "+  bestfit_table_name + " FROM '"+cwd+'/'+bestfitfile+"' DELIMITER ' ';"
    print(cmd)
    cursor.execute (cmd)
    db.commit()

    #Create the bestfit patch table
    cursor.execute ("DROP TABLE IF EXISTS "+bestfit_patch_table_name)
    print("DROP TABLE IF EXISTS "+bestfit_patch_table_name)
    cmd = "CREATE TABLE " + bestfit_patch_table_name + " (PatchID int, dmagcal double precision);"
    print(cmd)
    cursor.execute (cmd)
    db.commit()

    ##Copy the bestfit stars stuff in
    cmd = "COPY "+  bestfit_patch_table_name + " FROM '"+cwd+'/'\
          +bestfitpatchfile+"' DELIMITER ' ';"
    print(cmd)
    cursor.execute (cmd)
    db.commit()

    cursor.execute('DROP INDEX IF EXISTS starid_bf'+index_tag+' ;')
    #cursor.execute('DROP INDEX IF EXISTS master_patchid'+index_tag+' ;') #I think this is a typo with "master" there

    cmd = "CREATE INDEX patchid_patch" +index_tag +" on " + bestfit_patch_table_name + " (PatchID);"
    cursor.execute (cmd)

    cmd = "CREATE INDEX starid_bf" + index_tag +" on " +  bestfit_table_name + " (StarID);"
    cursor.execute (cmd)
    
    db.commit()



    
    print('finished, closing connection')
    cursor.close()
    db.close()



###################
#Easy way to load up a saved dictionary
def readDict(filename):
    d=n.load(filename)
    result={}
    for i in range(len(d.keys())):
        result[d.keys()[i]]=d[d.keys()[i]]
    return result

#Print the corected star_obs file
###################
def print_obsfile(ofile, starobs_output="star_obs_corrected.dat", stars=None):
    # This is the file that goes to selfCalib
    if ofile==None:
        ofile = open(starobs_output, 'w') 
        print("%s %s %s %s" %("#PatchID", "StarID", "StarObsMag", "StarMagObsErr"), file=ofile)
    if stars!=None:
        for obj in range(0, len(stars['id'])):
            print("%d %d %f %f" %(stars['fullpatch'][obj], stars['id'][obj],
                                           stars['rmagobs'][obj], stars['magerr'][obj]), file=ofile)
    return ofile

#Bin the stars spatially and by color to find delta mag terms
###################
def makeDeltamap(table1, nspace_bins=10, ncolor_bins=10, \
                 color_noise=0.05, delta_file='delta_data.npz', cc=False):
    """Generate a color-correction map using the  """

    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                               dbtype='pgsql', username=None, passwdname=None)

    #find the field of view max and mins
    cmd='select min(starx),min(stary),max(starx),max(stary) from '+table1+' ;'
    cursor.execute(cmd)
    field_sizes=cursor.fetchone()
    xsize=field_sizes[2]-field_sizes[0]
    ysize=field_sizes[3]-field_sizes[1]

    #set up the spatial bins
    xbins=(n.arange(nspace_bins+1)*xsize)/nspace_bins+field_sizes[0]
    ybins=(n.arange(nspace_bins+1)*ysize)/nspace_bins+field_sizes[1]

    #set up the color bins
    cmd = 'select min(StarTrueCol),max(StarTrueCol) from '+table1+' ;'
    cursor.execute(cmd)
    color_range=cursor.fetchone()
    color_bins= n.arange(ncolor_bins+1)* \
                (color_range[1]-color_range[0])/ncolor_bins+color_range[0]

    #arrays to hold the average delta_rs
    delta_rs=n.zeros([nspace_bins,nspace_bins,ncolor_bins])
    delta_rs_sigma=n.zeros([nspace_bins,nspace_bins,ncolor_bins])
    delta_rs_n=n.zeros([nspace_bins,nspace_bins,ncolor_bins])

    tag1=''
    tag2=''
    if cc == True:
        tag1='cc_obs'
        tag2='2'
    keys=['starid','StarObsMag','StarMagErr','StarTrueCol',\
          'dmagcal','starmagcalib']

    #loop over different radial bins (or x/y bins)
    for i in range(len(xbins)-1):
        for j in range(len(ybins)-1):
            cmd='select b.starid,b.StarObsMag,b.StarMagErr,\
            b.StarTrueCol, p.dmagcal, f.starmagcalib \
            from '+ table1 +' as b, '+ table1 + '_bestfitpatch'+tag2+' as p, \
            '+table1+'_bestfit'+tag2+' as f where b.patchid = p.patchid \
            and b.starid=f.starid and \
            b.starx >= %f and b.starx <= %f and b.stary >= %f and b.stary <= %f ;' \
            %(xbins[i], xbins[i+1], ybins[j],ybins[j+1])
        
            sqlresults=ui.sqlQuery(cursor,cmd, verbose=False)
            data=ui.assignResults(sqlresults,keys)
        
            if data['StarObsMag'] != None:
            
                data['delta_r']=data['StarObsMag']-data['dmagcal']-data['starmagcalib']
                #add noise to the color
                measured_color=data['StarTrueCol'] + \
                                color_noise*n.random.randn(len(data['StarTrueCol']))
                print("generating color map", i,j)
                for k in range(len(color_bins)-1):
                    in_bin=n.where((measured_color >= color_bins[k]) & \
                           (measured_color <= color_bins[k+1]) )
                    #feed the results into the arrays
                    delta_rs[i,j,k]=n.mean(data['delta_r'][in_bin])
                    delta_rs_sigma[i,j,k]=n.std(data['delta_r'][in_bin])
                    delta_rs_n[i,j,k]=len(data['delta_r'][in_bin])




    #set the center of the filter to be zero correction.
    x_center = n.where( n.abs(xbins) == n.min(n.abs(xbins)) )
    y_center = n.where( n.abs(ybins) == n.min(n.abs(ybins)) )
                    
    for k in range(len(color_bins)-1):
        delta_rs[:,:,k]=delta_rs[:,:,k]-max(delta_rs[x_center,y_center,k])


    n.savez(delta_file,delta_rs=delta_rs, delta_rs_sigma=delta_rs_sigma, \
                delta_rs_n=delta_rs_n, color_bins=color_bins, xbins=xbins,\
                ybins=ybins)
    
    #close the connection
    ui.sqlEndConnect(conn, cursor, verbose=True)

###################
def applyDeltamap(table1, delta_file='delta_data.npz',color_noise=0.05):

    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                               dbtype='pgsql', username=None, passwdname=None)

    keys=['fullpatch','id','rmagobs','magerr','true_color']
    #start the star_obs file
    ofile=print_obsfile(None, starobs_output='star_obs_corrected.dat')

    #load up the delta_table info

    delta_map=n.load(delta_file)
    xbins=delta_map['xbins']
    ybins=delta_map['ybins']
    color_bins=delta_map['color_bins']
    delta_rs=delta_map['delta_rs']

    for i in range(len(delta_map['xbins'])-1):
        for j in  range(len(delta_map['ybins'])-1):
            cmd='select b.patchid, b.starid,b.StarObsMag,b.StarMagErr,b.StarTrueCol \
                 from '+ table1 +' as b where  b.starx >= %f and b.starx <= %f and b.stary >= %f and b.stary <= %f ;' \
                 %(xbins[i], xbins[i+1], \
                   ybins[j],ybins[j+1])
            sqlresults=ui.sqlQuery(cursor,cmd, verbose=False)
            data=ui.assignResults(sqlresults,keys)
            if data['rmagobs'] != None:
                measured_color = data['true_color'] + \
                                 color_noise*n.random.randn(len(data['true_color']))
                #Note, not the same noise I added before, but probably fine
                
                measured_color.clip(min=n.min(color_bins), max=n.max(color_bins))
                del data['true_color']
                print('applying color map',i,j)

                for k in range(len(delta_map['color_bins'])-1):
                    in_bin=n.where((measured_color >= color_bins[k]) \
                                   & (measured_color <= color_bins[k+1]) )
                    data['rmagobs'][in_bin]=data['rmagobs'][in_bin]-\
                                             n.max(delta_rs[i,j,k])

                ofile=print_obsfile(ofile, stars=data)
    #close the connection
    ui.sqlEndConnect(conn, cursor, verbose=True)
    #close the ofile
    ofile.close()


def checkColor(stardata='stardata.dat',patchdata='patchdata.dat',\
               bestfit_star='test2_bestfit_Star.dat', title='Spatial Uniformity:', file=None):

    #read in the star files
    stardat=cp.read_true_stars(stardata)
    starcal=cp.read_bestfit_stars(bestfit_star)

    #match the stars up
    stardat = cp.match_stars_calvstrue(stardat, starcal)
    # adjust so they have the right overall zeropoint
    stardat = cp.adjust_stars_calvstrue(stardat)
    good=n.where(n.abs(stardat['magdiff']) < 0.05)
    if n.size(good) == 0:
        good=n.array([1])
    print('number of good stars = %f'%n.size(good))
    print('number of crazy bad fit stars = %f' %(n.size(stardat['magdiff']) - n.size(good)))
    print('RMS of stellar residuals = %f'%stardat['magdiff'][good].std())

    rms=stardat['magdiff'][good].std()

    pyl.figure()
    num, bins, patches = pyl.hist(stardat['magdiff'][good]*1000, bins=50)
    pyl.xlabel('m$_{true}$-m$_{best}$ (mmag)')
    pyl.title(title+', RMS=%.1f, clipped %.0f'%(rms*1000, n.size(stardat['magdiff']) - n.size(good)))
    pyl.savefig('StatsTruemBest.png', format='png')

    if file != None:
        data = stardat['magdiff'] +0
        #stats = {'n_stars':len(data), 'n_stars_used':len(good), 'max':n.max(data[good]), 'min':n.min(data[good]), 'mean':n.mean(data[good]),  'median':n.median(data[good]), 'std':n.std(data[good])}
        #ui.writeDatafile(file, stats, stats.keys())
        n_stars=len(data)
        n_stars_used=len(n.ravel(good))
        max=n.max(data[good])
        min=n.min(data[good])
        mean=n.mean(data[good])
        median=n.median(data[good])
        std=n.std(data[good])
        mean_all = n.mean(data)
        std_all = n.std(data)
        median_all = n.median(data)
        max_all= n.max(data)
        min_all = n.min(data)
        n.savez(file, n_stars=n_stars,n_stars_used=n_stars_used,max=max,min=min,mean=mean,median=median,std=std, \
                mean_all=mean_all, std_all=std_all, median_all=median_all, max_all=max_all, min_all=min_all)


def deltagrid(delta, cbin):

    nbins=n.size(delta['xbins'])
    scale=n.sqrt(n.max(delta['ybins'])**2+n.max(delta['xbins'])**2)
    x=(delta['xbins'][1:nbins]+delta['xbins'][0:nbins-1])/2./scale
    y=(delta['ybins'][1:nbins]+delta['ybins'][0:nbins-1])/2./scale
    ncb=n.size(delta['color_bins'])
    xx=x.repeat(10).reshape(10,10).transpose()
    yy=y.repeat(10).reshape(10,10)
    #goodi=n.where(delta['delta_rs_n'][:,:,cbin] < 10)
    #xx[goodi]=5000
    z=delta['delta_rs'][:,:,cbin]
    return xx,yy,z

def plotDelta(deltafile='delta_data.npz'):
    delta=readDict(deltafile)
    colors=[0,3,6,9]
    pyl.figure()
    pyl.subplots_adjust(wspace=0.3,hspace=0.3)
    for i in range(4):
        x,y,z=deltagrid(delta, n.max(colors[i]))
        #z=delta['delta_rs'][:,:,colors[i]]
        
        pyl.subplot(220+i+1)
        #pyl.contour(x,y,z,50)
        #pyl.savefig('delta.png', format='png')

        pyl.plot(n.sqrt(x**2+y**2),z*1000,'ko')
        pyl.xlabel('Radius (FoV')
        pyl.ylabel('$\Delta$r (mmag)')
        pyl.title('$g-i=$%f'%delta['color_bins'][colors[i]])

    pyl.savefig('delta.eps', format='eps')

    pyl.figure()
    x=n.ravel(x)
    y=n.ravel(y)
    z=n.ravel(z)
    ack=n.ravel(delta['delta_rs_n'][:,:,colors[0]])
    good=n.where(ack > 10)
    plot_contour_irregulargrid(x[good],y[good],z[good], nxbins=100, nybins=100)
    pyl.xlim([n.min(x),n.max(x)])
    pyl.ylim([n.min(y),n.max(y)])
    #pyl.contour(x,y,z,50)
    #pyl.title('$\Delta$r spatial map')
    pyl.savefig('delta_contour.eps', format='eps')


def find_zps(cursor,starid=5,patchid=5,patch_name='patch_name', \
             bestfit_name='bestfit_name'):
    """find the stupid magnitude zero points"""
    
    zj=n.zeros(len(starid))
    bf=n.zeros(len(starid))
    for i in range(len(zj)):
        cmd='select dmagcal from ' + patch_name + ' where PatchID = %i' %patchid[i]
        cursor.execute(cmd)
        zj[i]=n.max(cursor.fetchone())
        cmd = 'select  StarMagCalib from '+ \
          bestfit_name + ' where StarID = %i' %starid[i]
        cursor.execute(cmd)
        bf[i]=n.max(cursor.fetchone())

    return zj,bf

def plotStar(table, starid='4'):
    """Make a smimple plot of a single star to check how the sim is doing """

    name=table
    master_table_name=name
    patch_name = name+'_bestfitpatch'
    bestfit_name=name+'_bestfit'

    #open the connection the the DB
    #db = pg.connect(dsn = 'localhost:calsim')
    db = pg.connect(host='localhost', dbname='calsim')
    cursor = db.cursor ()
    cursor.execute('select max(starx),max(stary) from '+name+' ;')
    fieldscale=max(cursor.fetchone())


    cursor.execute('select (starx^2+stary^2)^.5,starid,startruecol,StarObsMag,PatchID from ' +name+' where starID='+starid)
    star1=n.array(cursor.fetchall () )

    star1_zps, star1_bf = find_zps(cursor, star1[:,1],star1[:,4],patch_name,bestfit_name)
    num,bins,patches = pyl.hist(star1[:,0]/fieldscale, histtype='bar', rwidth=0.8)
    pyl.xlabel('Radius from Center of FoV')
    pyl.ylabel('# of observations')

    pyl.savefig(name+'_hist.png',format='png')
    pyl.figure()

    pyl.plot(star1[:,0]/fieldscale,star1[:,3]-star1_zps-star1_bf,'ko')
    pyl.xlabel('Radius from Center of FoV' )
    pyl.ylabel('r$_{obs}$-zp-bestfit')
    pyl.title('Single star observations, g-i=%5f ' %star1[0,2])
    pyl.axhline(0.,linewidth=2, color='r')


    pyl.savefig(name+'singlestar.png',format='png')
    pyl.figure()


    #close the db connection
    print('finished, closing connection')
    cursor.close()
    db.close()


def iteratecc(table, index_tag):

    #put the updated fit into the database
    sim2db(table, index_tag, masterfile=None, \
           bestfitfile='test2_bestfit_Star.dat', bestfitpatchfile='test2_bestfit_Patch.dat')

    makeDeltamap(table)
    plotDelta()
    applyDeltamap(table)
    
#put the second round of bestfit parameters into the database
def corrected2db(name, index_tag, bfnum='2'):
    bestfit_table_name = name +'_bestfit'+bfnum
    bestfitfile='test'+bfnum+'_bestfit_Star.dat'
    bestfit_patch_table_name = name+'_bestfitpatch'+bfnum
    bestfitpatchfile='test'+bfnum+'_bestfit_Patch.dat'

    
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                               dbtype='pgsql', username=None, passwdname=None)

    cwd=os.getcwd()

    os.system("perl -i -p -e 's/^\s+//;' "+bestfitfile) ##remove leading spaces
    os.system("perl -i -p -e 'tr/ //s;' "+bestfitfile)  ##remove  any extra spaces
    os.system("perl -i -p -e 's/^\s+//;' "+bestfitpatchfile) ##remove leading spaces
    os.system("perl -i -p -e 'tr/ //s;' "+bestfitpatchfile)  ##remove  any extra spaces

    print("DROP TABLE IF EXISTS "+bestfit_table_name)    
    cursor.execute ("DROP TABLE IF EXISTS "+bestfit_table_name)
    cmd = "CREATE TABLE " + bestfit_table_name + " (StarID int, StarMagCalib double precision);"
    print(cmd)
    cursor.execute(cmd)

    cmd = "COPY "+  bestfit_table_name + " FROM '"+cwd+'/'+bestfitfile+"' DELIMITER ' ';"
    print(cmd)
    cursor.execute (cmd)

    #Create the bestfit patch table
    cursor.execute ("DROP TABLE IF EXISTS "+bestfit_patch_table_name)
    print("DROP TABLE IF EXISTS "+bestfit_patch_table_name)
    cmd = "CREATE TABLE " + bestfit_patch_table_name + " (PatchID int, dmagcal double precision);"
    print(cmd)
    cursor.execute (cmd)

    ##Copy the bestfit patch stuff in
    cmd = "COPY "+  bestfit_patch_table_name + " FROM '"+cwd+'/'\
          +bestfitpatchfile+"' DELIMITER ' ';"
    print(cmd)
    cursor.execute (cmd)

    cursor.execute('DROP INDEX IF EXISTS starid_bf'+index_tag+bfnum+' ;')
    cursor.execute('DROP INDEX IF EXISTS patchid_patch'+index_tag+bfnum+' ;')


    cmd = "CREATE INDEX patchid_patch" +index_tag +bfnum+" on " + bestfit_patch_table_name + " (PatchID);"
    cursor.execute (cmd)

    cmd = "CREATE INDEX starid_bf" + index_tag +bfnum+" on " +  bestfit_table_name + " (StarID);"
    cursor.execute (cmd)
 
    #close the connection
    ui.sqlEndConnect(conn, cursor, verbose=True)


#put the damn color corrected observed mags into the db
def color_corrected2db(name, num, filename='star_obs_corrected.dat'):
    cwd=os.getcwd()
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                               dbtype='pgsql', username=None, passwdname=None)

    #Create the bestfit stars table
    cursor.execute ("DROP TABLE IF EXISTS "+name+'cc_obs'+num)
    print("DROP TABLE IF EXISTS "+name+'cc_obs'+num)
    cmd = "CREATE TABLE " +  name+'cc_obs'+num + " (PatchID int, StarID int, StarObsMag double precision, StarMagObsErr double precision);"
    print(cmd)
    cursor.execute (cmd)

    #populate table
    cmd="COPY "+ name+'cc_obs'+num  + " FROM '"+cwd+'/'\
          +filename+"' DELIMITER ' ' HEADER CSV;"
    print(cmd)
    cursor.execute (cmd)

    #index things
    cmd = "CREATE INDEX starid_cc" + name+num+" on " + name+'cc_obs'+num  + " (StarID);"
    print(cmd)
    cursor.execute (cmd)

    #close the connection
    ui.sqlEndConnect(conn, cursor, verbose=True)
 


def bfpbsplot(name, bfnum='2', from_main=False, file=None):
    """Plot some histograms showing how well the fits are doing compared to reality """
    #grab from the db the star mag and bestfit patch zp
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                               dbtype='pgsql', username=None, passwdname=None)


    keys=['dr']

    t1=time.time()
    cc_obs='cc_obs'
    if (from_main == True):
        cc_obs=''

    #find out how many results to expect
    #print 'find the number of results to expect'
    cmd = 'select count(starid) from '+name+cc_obs
    cursor.execute(cmd)
    nrow=n.ravel(n.array(cursor.fetchall()))
    
    cmd='select m.starobsmag-p.dmagcal-b.starmagcalib FROM '+name+cc_obs+ \
         ' as m, '+name+'_bestfitpatch'+bfnum+' as p, '+name+'_bestfit'+bfnum+' as b WHERE m.patchid = p.patchid and m.starID = b.starID ;'

    cursor.execute(cmd)
    #print 'query complete, running fetchall', time.time()-t1

    #silly loop to make sure the memory usage doesn't explode when fetching the results from the DB
    if nrow < 1e6:
        data=cursor.fetchall()
        data = n.ravel(n.array(data))
    else:
        #decide how many fetch calls to make
        nf=int(n.floor(nrow/1e6)-1)
        d1=cursor.fetchmany(size = int(1e6)) #execute an initial call
        data=n.ravel(n.array(d1))   #start the data array
        for i in range (1,nf):
            #print i, time.time()-t1
            d1=cursor.fetchmany(size = int(1e6))
            data=n.concatenate((data, n.ravel(n.array(d1))))
        d1=cursor.fetchall() #fetch the remaining rows
        data=n.concatenate((data, n.ravel(n.array(d1))))
    good=n.where(abs(data) < 0.1)

    #now also grab the fitted patch minus true patch
    # read the data files from disk
    patchdat = cp.read_true_patchfile("patchdata.dat.s")
    patchcal = cp.read_bestfit_patchfile("test_bestfit_Patch.dat")

    # match the true and best-fit patch
    patchdat = cp.match_patch_calvstrue(patchdat, patchcal, sorted=True)

    # adjust the calibrated patch so they have approximately the right overall zeropoint
    patchdat = cp.adjust_patch_calvstrue(patchdat)
    ord=n.abs(patchdat['magdiff']).argsort()
    ord=ord[0:len(ord)*.95]
    patch_rms=n.std(patchdat['magdiff'][ord])
    goodpatch=n.where(abs(patchdat['magdiff']) < 0.1)
    #print 'generating histogram', time.time()-t1
    pyl.figure()
    ax = pyl.subplot(111)
    if n.size(good) > 1:
        num, bins, patches = pyl.hist(data[good]*1000,
                                  bins=50, normed=1, alpha=.75)
    else:
        bins = 50
    num, bins, patches = pyl.hist(patchdat['magdiff'][ord]*1000,
                                  bins=bins, normed=1, alpha=.75)
    #print 'adding labels etc', time.time()-t1
    pyl.xlabel('m$_{calibrated}$-m$_{best}$, patch$_{true}$-patch$_{best}$')
    rms=n.std(data[good])
    pyl.title(' RMS=%.1f mmag, patch RMS=%.1f'%(rms*1000,patch_rms*1000))
    #if there is a zero-patch, annotate with how it is offset
    nflux=0
    zeroPatch = n.where(patchdat['patchid'] == 0)
    zp_off=-666
    if n.size(zeroPatch) != 0:
        zp_off = patchdat['magdiff'][zeroPatch]
        print('The flux-star patch is %.2f mmag off from the best zero point'%(patchdat['magdiff'][zeroPatch]*1e3))
        pyl.text(0.02, 0.95,'$\Delta$ Flux ZP =%.2f mmag'%(patchdat['magdiff'][zeroPatch]*1e3), horizontalalignment='left',
                 verticalalignment='top', transform = ax.transAxes)
        cmd = 'select count(PatchID) from '+name+' where PatchID = 0;'
        cursor.execute(cmd)
        nflux = n.ravel(cursor.fetchall())
        nflux=nflux.astype('float')
        print('Number of Flux stars = %f'%nflux)
    pyl.savefig('StatsRepeat.png', format='png')

    ui.sqlEndConnect(conn, cursor, verbose=True)
    if file != None:
        good = n.where(data != -666)
        #stats = {'n_stars':len(data), 'n_stars_used':len(good), 'max':n.max(data[good]), 'min':n.min(data[good]), 'mean':n.mean(data[good]), 'median':n.median(data[good]), 'std':n.std(data[good])}
        #ui.writeDatafile(file, stats, stats.keys())
        n_stars=len(data)
        n_stars_used=len(n.ravel(good))
        max=n.max(data[good])
        min=n.min(data[good])
        mean=n.mean(data[good])
        median=n.median(data[good])
        std=n.std(data[good])
        patch_rms=patch_rms
        #n.savez(file, n_stars,n_stars_used,max,min,mean,median,std)
        n.savez(file, n_stars=n_stars,n_stars_used=n_stars_used,max=max,min=min,mean=mean,median=median,std=std, \
                patch_rms=patch_rms, nflux=nflux, zp_off=zp_off)

    

def star_sigma_clip(name, rms_clip=0.05, reweight=True):
    """Clip all observations of stars that have very large RMS values"""
    blocksize = 1e5 #largest number of stars to query for at a time
    ofile = print_obsfile(ofile=None, starobs_output = 'sig_clipped_star_obs.dat')
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                              dbtype='pgsql', username=None, passwdname=None)
    cmd = 'SELECT sd.starID, STDDEV(o.starObsMag - p.dmagcal - sd.startruemag) from '+name+'_bestfitpatch as p, '+name+' as o, '+name+'_stars as sd WHERE p.patchid = o.PatchID and o.starid = sd.starID GROUP BY sd.starid ORDER BY sd.starid;'
#    print cmd
#    cursor.execute(cmd)
    keys=['id','stds']
    data = ui.sqlQuery(cursor,cmd) #n.array(cursor.fetchall())
    data = ui.assignResults(data, keys)
    #print n.size(data)
    stds=data['stds']
    #starid=n.ravel(data[:,0])
    #med_std=n.median(stds)
    good=n.where( stds < rms_clip)#sig_clip*med_std)
    data['id'] = data['id'][good]
    nq = n.ceil((n.max(data['id'])+1)/blocksize)
    keys = ('fullpatch', 'id', 'rmagobs', 'magerr')
    for i in n.arange(nq):
        star_start=blocksize*i
        star_end = blocksize*(i+1)
        if star_start < n.max(data['id']):            
            cmd = 'SELECT PatchID, starID, starObsMag, StarMagErr from '+name+' WHERE starid >= %f and starid < %f'%(star_start, star_end)
            stars = ui.sqlQuery(cursor,cmd)
            stars = ui.assignResults(stars, keys)
            #now to match stars['id'] to starid and only keep the good ones
            #print 'data returned from db ', type(stars['id']), n.size(stars['id'])
            #print 'data where stars are good ', type(starid), n.size(starid)
            condition = n.in1d(stars['id'], data['id'])
            for key in keys:
                stars[key]=stars[key][condition]
            ofile = print_obsfile(ofile,stars = stars)
    ui.sqlEndConnect(conn, cursor, verbose=True)


def clipReweight(name, rms_clip=0.04, reweight=True, error_min=0.003):
    """Clip stars with high RMS and calc the RMS per patch and change the magerr to the patch RMS"""
    blocksize = 1e5 #number of patches to grab at once
    ofile = print_obsfile(ofile=None, starobs_output = 'sig_clipped_star_obs.dat')
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                              dbtype='pgsql', username=None, passwdname=None)
    #find the RMS per star, sorted by star ID
    cmd = 'SELECT sbf.starID, STDDEV(o.starObsMag - p.dmagcal - sbf.starmagcalib) from '+name+'_bestfitpatch as p, '+name+' as o, '+name+'_bestfit as sbf WHERE p.patchid = o.PatchID and o.starid = sbf.starID GROUP BY sbf.starid ORDER BY sbf.starid;'
    keys=['id','stds']
    print("Finding the RMS per star")
    stardata = ui.sqlQuery(cursor,cmd) #n.array(cursor.fetchall())
    stardata = ui.assignResults(stardata, keys)
    #clip off the high RMS stars
    condition = n.where(stardata['stds'] < rms_clip)
    print('clipping %i stars due to high RMS in repeat observations'%(n.size(stardata['stds'])-n.size(condition)))
    for key in stardata.keys():
        stardata[key]=stardata[key][condition]
    print('got stats for %i stars'%n.size(stardata['id']))
    #find the RMS per patch
    if reweight:        
        print("Finding RMS per patch")
        cmd = 'SELECT o.PatchID, STDDEV(o.starObsMag  - sd.startruemag) from '+name+'_bestfitpatch as p, '+name+' as o, '+name+'_stars as sd, '+name+'_bestfit as sbf WHERE p.patchid = o.PatchID and o.starid = sd.starID and o.starid = sbf.starid GROUP BY o.PatchID;'
        keys=['id','stds']
        patchdata = ui.sqlQuery(cursor,cmd)
        patchdata=ui.assignResults(patchdata,keys)
        low_error = n.where(patchdata['stds'] < error_min)
        patchdata['stds'][low_error]=error_min
        print('got stats on %i patches'%n.size(patchdata['stds']))
    #now to query for blocks of patchid, starid, starobsmag, starobsmagerr.  I think grab blocks of patches
    cmd='SELECT max(patchid),min(patchid) from '+name
    cursor.execute(cmd)
    patchRange = n.ravel(cursor.fetchall())
    num = n.ceil((n.max(patchRange)-n.min(patchRange))/blocksize)
    keys = ('fullpatch', 'id', 'rmagobs', 'magerr')
    print("grabbing %i chunks of patches and reweighting errors by patch RMS"%num)
    for i in n.arange(num):
        minp = n.min(patchRange)+blocksize*i
        maxp = n.min(patchRange)+blocksize*(i+1)
        cmd = 'SELECT PatchID, starID, starObsMag, StarMagErr from '+name+' WHERE PatchID >= %f and PatchID < %f order by PatchID'%(minp, maxp)
        stars = ui.sqlQuery(cursor,cmd)
        stars = ui.assignResults(stars,keys)
        condition = n.in1d(stars['id'], stardata['id'])
        for key in keys:
            stars[key]=stars[key][condition]
        if reweight:
            left = n.searchsorted(stars['fullpatch'],patchdata['id'], side='left')
            right = n.searchsorted(stars['fullpatch'],patchdata['id'], side='right')
            for j in n.arange(n.size(left)):
                stars['magerr'][left[j]:right[j]] = patchdata['stds'][j]
        ofile = print_obsfile(ofile, stars = stars)
    ofile.close()
    ui.sqlEndConnect(conn, cursor, verbose=True)

def star_out_clip(name, outclip=0.05, dbname='calsim', reweight=True):
    """Clip off the major outliers for individual observations """
    blocksize = 1e6 #largest number of stars to query for at a time
    ofile = print_obsfile(ofile=None, starobs_output = 'out_clipped_star_obs.dat')
    keys = ('fullpatch', 'id', 'rmagobs', 'magerr', 'resid')
    conn,cursor=ui.sqlConnect('localhost', dbname=dbname,\
                              dbtype='pgsql', username=None, passwdname=None)

    cmd = 'select count(starid) from '+name
    cursor.execute(cmd)
    nrows=n.ravel(cursor.fetchall())
    cmd='select m.patchid, m.starid, m.starobsmag, m.starmagerr, m.starobsmag-p.dmagcal-b.starmagcalib FROM '+name + \
         ' as m, '+name+'_bestfitpatch as p, '+name+'_bestfit as b WHERE m.patchid = p.patchid and m.starID = b.starID order by m.patchid;'
    results = massivedb(cursor, cmd,nrows,keys)
    good = n.ravel(n.where(n.abs(results['resid']) < outclip))
    bad = n.ravel(n.where(n.abs(results['resid']) > outclip))
    #need to clobber the old star data file, and delete the right rows from master table so the later analysis works.
    print('deleting %i observations from the master database, leaving %i stars'%(n.size(bad),n.size(good)))
    for i in n.arange(n.size(bad)):
        cmd = 'delete from '+name+' where starid = %d and patchid = %d'%(results['id'][bad[i]],results['fullpatch'][bad[i]])
        #print i,cmd
        cursor.execute(cmd)
    for key in results.keys():
        results[key] = results[key][good]
    print('resulting obsfile should have %i stars in it'%(n.size(results['id'])))
    if reweight:
        cmd = 'SELECT o.PatchID, STDDEV(o.starObsMag - p.dmagcal - sd.startruemag) from '+name+'_bestfitpatch as p, '+name+' as o, '+name+'_stars as sd, '+name+'_bestfit as sbf WHERE p.patchid = o.PatchID and o.starid = sd.starID and o.starid = sbf.starid GROUP BY o.PatchID order by o.patchid;'
        print(cmd)
        cursor.execute(cmd)
        data=n.array(cursor.fetchall())
        patch_id = n.ravel(data[:,0]) 
        patch_resid_std = n.ravel(data[:,1])
        upatch = n.unique(results['fullpatch'])
        left = n.searchsorted(results['fullpatch'],patch_id, side='left')
        right=n.searchsorted(results['fullpatch'],patch_id, side='right')
        for i in n.arange(n.size(left)):
            results['magerr'][left[i]:right[i]]=patch_resid_std[i]
        #for patch in upatch:
        #    new_err = patch_resid_std[n.where(patch_id == patch)]
        #    if n.max(new_err) != None:
        #        match_patch = n.ravel(n.where(results['fullpatch']==patch))
        #        results['magerr'][match_patch]=match_patch*0.+n.max(new_err)
    ofile = print_obsfile(ofile, stars = results)
    ofile.close()
    ui.sqlEndConnect(conn, cursor, verbose=True)

def plotperstar(name, file=None):
    """Plot up some histograms where the values have been averaged for each starID  """
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                              dbtype='pgsql', username=None, passwdname=None)
#    cmd = 'SELECT sd.starID, AVG(o.starObsMag - p.dmagcal - sd.startruemag), STDDEV(o.starObsMag - p.dmagcal - sd.startruemag) from '+name+'_bestfitpatch as p, '+name+' as o, '+name+'_stars as sd WHERE p.patchid = o.PatchID and o.starid = sd.starID GROUP BY sd.starid;'
    cmd = 'SELECT sd.starID, AVG(o.starObsMag - p.dmagcal - sd.startruemag), STDDEV(o.starObsMag - p.dmagcal - sd.startruemag) from '+name+'_bestfitpatch as p, '+name+' as o, '+name+'_stars as sd, '+name+'_bestfit as sbf WHERE p.patchid = o.PatchID and o.starid = sd.starID and o.starid = sbf.starid GROUP BY sd.starid;'
    print(cmd)
    cursor.execute(cmd)
    data = n.array(cursor.fetchall())
    floating_zp = n.median(data[:,1])
    data[:,1]=data[:,1]-floating_zp
    nstars=len(data[:,1])
    good=n.where(n.abs(data[:,1]) < 0.1)
    avs=n.ravel(data[good,1])
    stds=n.ravel(data[:,2])
    stds=stds[n.where((stds > 0) & (stds < .5))]
    del data

    pyl.figure()
    num,bins,patches = pyl.hist(avs*1000,bins=50)
    pyl.xlabel('<m$_{calibrated}$-m$_{true}$>$_i$ (mmag)')
    pyl.title('RMS=%.1f mmag, clipped %i stars'
              %(n.std(avs*1000),nstars-len(avs)) )
    pyl.savefig('StatsAveStar.png',format='png')

    pyl.figure()
    num,bins,patches = pyl.hist(stds*1000,bins=50)
    if n.max(bins) > 150:
        pyl.figure()
        num,bins,patches = pyl.hist(stds*1000,bins=500, log=True)       
    pyl.xlabel('RMS(m$_{calibrated}$-m$_{true}$)$_i$ (mmag)')
    pyl.title('Repeatability: Median=%.1f mmag, clipped %i stars'
              %(n.median(stds*1000), nstars-len(stds)))

    pyl.savefig('StatsRMSStars.png',format='png')

    #make a cummultative plot for the RMS values
    stds.sort()
    y=n.arange(len(stds),dtype='float')
    y=y/n.max(y)*100.
    fifty_level=stds[n.min(n.where(y > 50))]
    ninety_level=stds[n.min(n.where(y > 90))]
    p='%'
    s1='50'+p+' level = %.1f, 90'%(fifty_level*1000)+p
    s2=' level = %.1f'%(ninety_level*1000)
    pyl.figure()
    pyl.plot(stds*1000,y)
    #make some lines at 50% and 90%
    pyl.xlabel('RMS(m$_{calibrated}$-m$_{true}$)$_i$ mmag')
    pyl.ylabel('cumulative percentage')
    pyl.axhline(y=50, linestyle='--')
    pyl.axhline(y=90, linestyle='--')
   
    pyl.title(s1+s2)
    pyl.savefig('StatsCumRMS.png',format='png')

    #grab the stats per patch
    cmd = 'SELECT o.PatchID, AVG(o.starObsMag - p.dmagcal - sd.startruemag), STDDEV(o.starObsMag - p.dmagcal - sd.startruemag) from '+name+'_bestfitpatch as p, '+name+' as o, '+name+'_stars as sd, '+name+'_bestfit as sbf WHERE p.patchid = o.PatchID and o.starid = sd.starID and o.starid = sbf.starid GROUP BY o.PatchID;'
    print(cmd)
    cursor.execute(cmd)
    data=n.array(cursor.fetchall())
    patch_resid_ave = n.ravel(data[:,1]- n.ravel(floating_zp)) #remove floating zp
    patch_resid_std = n.ravel(data[:,2])
    print('a few elements from the std array', patch_resid_std[0:4])
    del data
    pyl.figure()
    num,bins,patches = pyl.hist(patch_resid_ave[n.where(n.abs(patch_resid_ave) < 0.05 )]*1000, bins=50)
    pyl.xlabel('<m$_{calibrated}$-m$_{true}$>$_j$ mmag')
    pyl.ylabel('#')
    pyl.title('Error Averaged Per Patch, RMS=%.0f'%(n.std(patch_resid_ave[n.where(n.abs(patch_resid_ave) < 0.05 )])*1000))
    pyl.savefig('StatsPatchRMS.png', format='png')
    pyl.figure()
    good = n.where(patch_resid_std > 0.)
    num,bins,patches = pyl.hist(patch_resid_std[good]*1000, bins=50)  #[n.where(n.abs(patch_resid_ave) < 0.05 )]
    pyl.xlabel('$\sigma_j^{patch}$ mmag')
    pyl.ylabel('#')
    nclip = n.size(patch_resid_ave)-n.size(n.where(n.abs(patch_resid_ave) < 0.05 ))
    pyl.title('RMS Per Patch, median=%.0f mmag'%(n.median(patch_resid_std[good]*1000))) #, clipped = %i'%(n.median(patch_resid_std[n.where(n.abs(patch_resid_ave) < 0.05 )])*1000),nclip)
    pyl.savefig('StatsPatchRMS2.png', format='png')

    ui.sqlEndConnect(conn, cursor, verbose=True)

    if file != None:
        data = stds+0
        good = n.where(data != -666)
        #stats = {'n_stars':len(data), 'n_stars_used':len(n.ravel(good)), 'max':n.max(data[good]), 'min':n.min(data[good]), 'mean':n.mean(data[good]), 'median':n.median(data[good]), 'std':n.std(data[good])}
        n_stars=len(data)
        n_stars_used=len(n.ravel(good))
        max=n.max(data[good])
        min=n.min(data[good])
        mean=n.mean(data[good])
        median=n.median(data[good])
        std=n.std(data[good])
        patch_std = n.median(patch_resid_std)
        patch_ave = n.median(patch_resid_ave)
        patch_ave_ssqq = n.median( (patch_resid_ave*patch_resid_ave)**0.5)
        patch_ave_spread = n.std(patch_resid_ave[n.where(n.abs(patch_resid_ave) < 0.05 )])
        levels=( ('50', fifty_level), ('90',ninety_level))
        #n.savez(file, n_stars,n_stars_used,max,min,mean,median,std)
        n.savez(file, n_stars=n_stars,n_stars_used=n_stars_used,\
                max=max,min=min,mean=mean,median=median,std=std, levels=levels, patch_std=patch_std, patch_ave=patch_ave,patch_ave_ssqq=patch_ave_ssqq, patch_ave_spread=patch_ave_spread)
        #ui.writeDatafile(file, stats, stats.keys())

def plotResidmap(name,nx=10, ny=10):
    """Use the database to divide up the focal plane and look at the residuals"""
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                               dbtype='pgsql', username=None, passwdname=None)

    #find the field of view max and mins
    cmd='select min(starx),min(stary),max(starx),max(stary) from '+name+' ;'
    cursor.execute(cmd)
    field_sizes=cursor.fetchone()
    xsize=field_sizes[2]-field_sizes[0]
    ysize=field_sizes[3]-field_sizes[1]
    
    #set up the spatial bins
    xbins=(n.arange(nx+1)*xsize)/nx+field_sizes[0]
    ybins=(n.arange(ny+1)*ysize)/ny+field_sizes[1]
    
    #find the maximum radius that the stars are at
    #cmd='select max(sqrt(starx^2+stary^2)) from '+name+' ;'
    #cursor.execute(cmd)
    #max_rad=cursor.fetchone()
    max_rad=0.03142675

    #generate dictionary to hold results with number of stars, mean, median, stdev
    result={'xmin':n.zeros([nx,ny]),'xmax':n.zeros([nx,ny]),'ymin':n.zeros([nx,ny]),\
            'ymax':n.zeros([nx,ny]), 'nstars':n.zeros([nx,ny]), 'mean':n.zeros([nx,ny]),\
            'median':n.zeros([nx,ny]),'std':n.zeros([nx,ny])}

    #loop over the spatial bins
    #want the observed mag - any color correction - patch zp - bestfit mag
    for i in range(len(xbins)-1):
        for j in range(len(ybins)-1):
            if (xbins[i]**2+ybins[j]**2 < max_rad**2) & (xbins[i+1]**2+ybins[j+1]**2 < max_rad**2)\
                   & (xbins[i+1]**2+ybins[j]**2 < max_rad**2) & (xbins[i]**2+ybins[j+1]**2 < max_rad**2):
               print(i,j)
               cmd='select m.starobsmag-p.dmagcal-m.startruemag FROM '+name+\
                    ' as m, '+name+'_bestfitpatch as p '+\
                    'WHERE m.patchid=p.patchid  and m.starx > %f'%(xbins[i])+\
                    ' and m.starx < %f and m.stary > %f'%(xbins[i+1],ybins[j])+\
                    ' and m.stary < %f limit 10000 ;'%(ybins[j+1])
               cursor.execute(cmd)
               #now I think this should return a small enough number of results that I don't need to 
               temp_result=cursor.fetchall()
               result['xmin'][i,j]=xbins[i]
               result['xmax'][i,j]=xbins[i+1]
               result['ymin'][i,j]=ybins[j]
               result['ymax'][i,j]=ybins[j+1]
               result['nstars'][i,j]=len(temp_result)
               result['mean'][i,j]=n.mean(temp_result)
               result['median'][i,j]=n.median(temp_result)
               result['std'][i,j]=n.std(temp_result)


    #quick cludge to remove the overall zeropoint offset
    zp= n.median(result['median'][n.where(result['nstars'] != 0)])
    result['mean']=result['mean']-zp
    result['median']=result['median']-zp
    
    x=(result['xmin']+result['xmax'])/2
    y=(result['ymin']+result['ymax'])/2
    r=(x**2+y**2)**0.5
    r=r/n.max(r)
    #CS=pyl.contour(x,y,result['nstars'])
    #pyl.colorbar()
    #pyl.title('Number of stars per bin')
    #pyl.savefig('rmap_nstars.png',format='png')
    #pyl.figure()

    good=n.where(result['nstars'] != 0)
    #make some radial plots
    fig=pyl.figure()
    fig.subplots_adjust(wspace=0.4)

    pyl.subplot(221)
    pyl.plot(r[good],result['nstars'][good],'ko')
    pyl.xlabel('Radius')
    pyl.ylabel('# of stars per bin')

    pyl.subplot(222)
    pyl.plot(r[good],result['mean'][good]*1000., 'ko')
    pyl.xlabel('Radius')
    pyl.ylabel('Mean residual value (mmag)')


    pyl.subplot(223)
    pyl.plot(r[good],result['median'][good]*1000, 'ko')
    pyl.xlabel('Radius')
    pyl.ylabel('Median residual value (mmag)')

    pyl.subplot(224)
    pyl.plot(r[good],result['std'][good]*1000, 'ko')
    pyl.xlabel('Radius')
    pyl.ylabel('RMS residual value (mmag)')

    pyl.savefig('rmap_radplots.png',format='png')

    pyl.figure()
    pyl.imshow(result['std'], extent=[n.min(x),n.max(x), n.min(y),n.max(y)],cmap=pyl.cm.jet,\
               vmin=n.min(result['std'][n.where(result['nstars'] != 0)]))
    pyl.colorbar()
    pyl.savefig('rmap_space.png',format='png')

    ui.sqlEndConnect(conn, cursor, verbose=True)
    return result


def allStats(name, file):
    """grab some stats on the mag offsets from the database  """
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                              dbtype='pgsql', username=None, passwdname=None)
    cmd = 'SELECT  AVG(dMag_var),AVG(dMag_snr),AVG(dMag_zp),AVG(dMag_zp_dist),AVG(dMag_color),AVG(dMag_color_dist),AVG(dMag_sin),AVG(dMag_flat_sin),AVG(dMag_cloud),AVG(dMag_cloud_image),AVG(dMag_rr), AVG(dMag_kep), AVG(dMag_gainvar), AVG(dMag_exptvar) from '+name
    cursor.execute(cmd)
    avgs = n.array(cursor.fetchall())
    cmd = 'select STDDEV(dMag_var),STDDEV(dMag_snr),STDDEV(dMag_zp),STDDEV(dMag_zp_dist),STDDEV(dMag_color),STDDEV(dMag_color_dist),STDDEV(dMag_sin),STDDEV(dMag_flat_sin),STDDEV(dMag_cloud), STDDEV(dMag_cloud_image), STDDEV(dMag_rr), STDDEV(dMag_kep), STDDEV(dMag_gainvar), STDDEV(dMag_exptvar) from '+name
    cursor.execute(cmd)
    stds =  n.array(cursor.fetchall())
    cmd =  'select MIN(dMag_var),MIN(dMag_snr),MIN(dMag_zp),MIN(dMag_zp_dist),MIN(dMag_color),MIN(dMag_color_dist),MIN(dMag_sin),MIN(dMag_flat_sin),MIN(dMag_cloud),MIN(dMag_cloud_image),MIN(dMag_rr), MIN(dMag_kep), MIN(dMag_gainvar), MIN(dMag_exptvar) from '+name
    cursor.execute(cmd)
    mins= n.array(cursor.fetchall())
    
    cmd =  'select MAX(dMag_var),MAX(dMag_snr),MAX(dMag_zp),MAX(dMag_zp_dist),MAX(dMag_color),MAX(dMag_color_dist),MAX(dMag_sin),MAX(dMag_flat_sin),MAX(dMag_cloud),MAX(dMag_cloud_image),MAX(dMag_rr), MAX(dMag_kep), MAX(dMag_gainvar), MAX(dMag_exptvar) from '+name
    cursor.execute(cmd)
    maxes= n.array(cursor.fetchall())

    names = ['dMag_var','dMag_snr','dMag_zp','dMag_zp_dist','dMag_color','dMag_color_dist','dMag_sin',
             'dMag_flat_sin','dMag_cloud','dMag_cloud_image','dMag_rr', 'dMag_kep', 'dmag_gainvar','dmag_exptvar' ]

    n.savez(file, avgs=avgs, stds=stds,mins=mins,maxes=maxes,names=names)
    ui.sqlEndConnect(conn, cursor, verbose=True)

def flat_correct_pernight(sim_name, nx=10,ny=10, plotnights=[31]):
    #build a spatial flat field map for each night and apply it to generate a new star_obs.dat file
    ofile=print_obsfile(None, starobs_output='star_obs_flat.dat')
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                              dbtype='pgsql', username=None, passwdname=None)
    #get the number of nights that need to be looped over
    cmd='select distinct night from %s'%(sim_name)
    #cursor.execute(cmd)
    nights=ui.sqlQuery(cursor,cmd, verbose=False)
    nights=ui.assignResults(nights,['nights'])
    nights=nights['nights'] #that's a roundabout way to do that

    keys=['id', 'fullpatch','starx','stary','rmagobs', 'magerr','dm']
    for i in nights:
        #print n.where(nights == i)
        print(i)
        cmd='select a.starid,a.patchid,a.starx,a.stary,a.starobsmag,a.starmagerr,a.starobsmag-bfp.dmagcal-bf.starmagcalib from %s as a, %s_bestfit as bf, %s_bestfitpatch as bfp where bf.starid=a.starid and bfp.patchid=a.patchid and  a.night=%i;'%(sim_name,sim_name,sim_name,n.max(i))
        result=ui.sqlQuery(cursor,cmd, verbose=False)
        data_in=ui.assignResults(result,keys)
        #corrected_obsmag=n.copy(data_in['rmagobs'])
        x=y=z=num=n.array(0.)
        if data_in['starx'] != None:
            xbins=n.arange(nx+1)*(n.max(data_in['starx'])-\
                              n.min(data_in['starx']))/nx+\
                              n.min(data_in['starx'])
            ybins=n.arange(ny+1)*(n.max(data_in['stary'])-\
                              n.min(data_in['stary']))/ny+\
                              n.min(data_in['stary'])
            #make a map of the correction and apply it
            for j in range(nx):
                for k in range(ny):
                    inbin=n.where( (data_in['starx'] >= xbins[j]) \
                               & (data_in['starx'] <= xbins[j+1]) \
                               & (data_in['stary'] >= ybins[k]) \
                               & (data_in['stary'] <= ybins[k+1]))
                    if len(n.ravel(inbin)) > 3:  
                        #print 'dm',n.mean(data_in['dm'][inbin])
                        m=n.mean(data_in['dm'][inbin])
                        if n.abs(m)-3.*n.abs(m)/n.sqrt(len(n.ravel(inbin))) > 0.:
                             data_in['rmagobs'][inbin]=data_in['rmagobs'][inbin]\
                             -n.mean(data_in['dm'][inbin])
                             x=n.append(x,xbins[j])
                             y=n.append(y,ybins[k])
                             z=n.append(z,m)
                             num=n.append(num,n.size(inbin))
                             
            if i in plotnights:
                x=x[1:]
                y=y[1:]
                z=z[1:]
                num=num[1:]
                pyl.figure()
                plot_contour_irregulargrid(x,y,z, nxbins=10, nybins=10)
                pyl.title('Residual Correction Map')
                pyl.savefig('MapCorrect.png',format='png')
                pyl.clim(vmin=-4.5,vmax=4.5)
                pyl.figure()
                plot_contour_irregulargrid(x,y,num, nxbins=10, nybins=10)
                pyl.title('Number of Stars in Map')
                pyl.savefig('MapCorrectN.png',format='png')
        
            ofile=print_obsfile(ofile,stars=data_in)
    ui.sqlEndConnect(conn, cursor, verbose=True)

def flat_correct_nightblock(sim_name, nx=10,ny=10, night_blocks=30, plotblock=[0]):
    #build a spatial flat field map for blocks of nights and apply it to generate a new star_obs.dat file
    ofile=print_obsfile(None, starobs_output='star_obs_flat.dat')
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                              dbtype='pgsql', username=None, passwdname=None)
    #get the number of nights that need to be looped over
    cmd='select distinct night from %s'%(sim_name)
    #cursor.execute(cmd)
    nights=ui.sqlQuery(cursor,cmd, verbose=False)
    nights=ui.assignResults(nights,['nights'])
    nights=n.sort(nights['nights']) #that's a roundabout way to do that
    night_chunks=n.arange((n.max(nights)-n.min(nights))/night_blocks)*night_blocks
    if n.max(night_chunks) < n.max(nights):
        night_chunks=n.append(night_chunks,n.max(nights)+1)
    keys=['id', 'fullpatch','starx','stary','rmagobs', 'magerr','dm']
    for i in n.arange(n.size(night_chunks)-1):
        #print n.where(nights == i)
        print(i)
        night_min=night_chunks[i]+1 #this should works since night begins at 1.
        night_max=night_chunks[i+1]
        cmd='select a.starid,a.patchid,a.starx,a.stary,a.starobsmag,a.starmagerr,a.starobsmag-bfp.dmagcal-bf.starmagcalib from %s as a, %s_bestfit as bf, %s_bestfitpatch as bfp where bf.starid=a.starid and bfp.patchid=a.patchid and  a.night >= %f and a.night < %f;'%(sim_name,sim_name,sim_name,night_min,night_max)
        result=ui.sqlQuery(cursor,cmd, verbose=False)
        data_in=ui.assignResults(result,keys)
        #corrected_obsmag=n.copy(data_in['rmagobs'])
        x=y=z=num=n.array(0.)
        if data_in['starx'] != None:
            xbins=n.arange(nx+1)*(n.max(data_in['starx'])-\
                              n.min(data_in['starx']))/nx+\
                              n.min(data_in['starx'])
            ybins=n.arange(ny+1)*(n.max(data_in['stary'])-\
                              n.min(data_in['stary']))/ny+\
                              n.min(data_in['stary'])
            #make a map of the correction and apply it
            for j in range(nx):
                for k in range(ny):
                    inbin=n.where( (data_in['starx'] >= xbins[j]) \
                               & (data_in['starx'] <= xbins[j+1]) \
                               & (data_in['stary'] >= ybins[k]) \
                               & (data_in['stary'] <= ybins[k+1]))
                    if len(n.ravel(inbin)) > 3:  
                        #print 'dm',n.mean(data_in['dm'][inbin])
                        m=n.mean(data_in['dm'][inbin])
                        if n.abs(m)-3.*n.abs(m)/n.sqrt(len(n.ravel(inbin))) > 0.:
                             data_in['rmagobs'][inbin]=data_in['rmagobs'][inbin]\
                             -n.mean(data_in['dm'][inbin])
                             x=n.append(x,xbins[j])
                             y=n.append(y,ybins[k])
                             z=n.append(z,m)
                             num=n.append(num,n.size(inbin))
                             
            if i in plotblock:
                x=x[1:]
                y=y[1:]
                z=z[1:]
                num=num[1:]
                pyl.figure()
                plot_contour_irregulargrid(x,y,z*1000, nxbins=nx, nybins=ny)
                pyl.title('Residual Correction Map')
                pyl.savefig('MapCorrect%i.png'%i,format='png')
                pyl.figure()
                plot_contour_irregulargrid(x,y,num, nxbins=nx, nybins=ny)
                pyl.title('Number of Stars in Map')
                pyl.savefig('MapCorrectN%i.png'%i,format='png')
        
            ofile=print_obsfile(ofile,stars=data_in)
    ui.sqlEndConnect(conn, cursor, verbose=True)


def massivedb(cursor, cmd,nrows,keys):
    #retrieve results from sql querry in blocks so memory usage doesn't explode
    print('executing large sql')
    cursor.execute(cmd)
    nf=int(n.floor(nrows)/1e6)-1
    print('fetching first million')
    d1=n.array(cursor.fetchmany(size=int(1e6)))
    result={}
    for i in range(n.size(keys)):
        result[keys[i]]=d1[:,i]
    for i in range(1,nf):
        d1=n.array(cursor.fetchmany(size=int(1e6)))
        for j in range(n.size(keys)):
            result[keys[j]]=n.concatenate( (result[keys[j]],d1[:,j]))
    return result
    

def flat_correct_global(sim_name, nx=10,ny=10):
    #build a spatial flat field map for each night and apply it to generate a new star_obs.dat file
    ofile=print_obsfile(None, starobs_output='star_obs_flat.dat')
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                              dbtype='pgsql', username=None, passwdname=None)

    keys=['id', 'fullpatch','starx','stary','rmagobs', 'magerr','dm']
    cmd='select count(starx) from %s'%(sim_name)
    cursor.execute(cmd)
    count=n.ravel(cursor.fetchall())
    print('retrieved count')
    
    cmd='select a.starid,a.patchid,a.starx,a.stary,a.starobsmag,a.starmagerr,a.starobsmag-bfp.dmagcal-bf.starmagcalib from %s as a, %s_bestfit as bf, %s_bestfitpatch as bfp where bf.starid=a.starid and bfp.patchid=a.patchid ;'%(sim_name,sim_name,sim_name)

    data_in=massivedb(cursor,cmd,count,keys)
    x=y=z=num=n.array(0.)
    if data_in['starx'] != None:
        xbins=n.arange(nx+1)*(n.max(data_in['starx'])-\
                              n.min(data_in['starx']))/nx+\
                              n.min(data_in['starx'])
        ybins=n.arange(ny+1)*(n.max(data_in['stary'])-\
                              n.min(data_in['stary']))/ny+\
                              n.min(data_in['stary'])
        for j in range(nx):
            for k in range(ny):
                inbin=n.where( (data_in['starx'] >= xbins[j]) \
                               & (data_in['starx'] <= xbins[j+1]) \
                               & (data_in['stary'] >= ybins[k]) \
                               & (data_in['stary'] <= ybins[k+1]))
                if len(n.ravel(inbin)) > 3:  
                    m=n.mean(data_in['dm'][inbin])
                    if n.abs(m)-3.*n.abs(m)/n.sqrt(len(n.ravel(inbin))) > 0.:
                        data_in['rmagobs'][inbin]=data_in['rmagobs'][inbin]\
                             -n.mean(data_in['dm'][inbin])
                        x=n.append(x,xbins[j])
                        y=n.append(y,ybins[k])
                        z=n.append(z,m)
                        num=n.append(num,n.size(inbin))
        ofile=print_obsfile(ofile,stars=data_in)
    ui.sqlEndConnect(conn, cursor, verbose=True)
    x=x[1:]
    y=y[1:]
    z=z[1:]
    num=num[1:]
    pyl.figure()
    plot_contour_irregulargrid(x,y,z*1000, nxbins=10, nybins=10)
    pyl.clim(vmin=-4.5,vmax=4.5)
    pyl.title('Residual Correction Map')
    pyl.savefig('MapCorrect.png',format='png')
    pyl.figure()
    plot_contour_irregulargrid(x,y,num, nxbins=10, nybins=10)
    pyl.title('Number of Stars in Map')
    pyl.savefig('MapCorrectN.png',format='png')
      


def residCheck(table,night=None, visitID= None, dmag='dmag_var+dmag_snr+dmag_zp_dist+dmag_color+dmag_color_dist+dmag_sin+dmag_flat_sin+dmag_cloud+dmag_rr', fname='residmap.png'):
    """Plot up the applied residuals from one night """
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                               dbtype='pgsql', username=None, passwdname=None)

    if night == None:
        cmd = 'select min(night) from %s where night > 30'%table
        cursor.execute(cmd)
        night=n.array(cursor.fetchall(), dtype='int')
    
    if visitID != None:
        cmd='select distinct visitid from %s where night=%i limit 10'%(table,night)
        cursor.execute(cmd)
        uvids=n.array(cursor.fetchall(),dtype='int')
        print(uvids[0])
        cmd='select  starx,stary,%s from %s where night=%i and visitid=%i'%(dmag,table,night,uvids[0])
    else:
        cmd='select  starx,stary,%s from %s where night=%i'%(dmag,table,night)
    cursor.execute(cmd)
    keys=['starx','stary','dmag']
    sqlresults=ui.sqlQuery(cursor,cmd, verbose=False)

    data=ui.assignResults(sqlresults,keys)

    #pyl.plot(data['starx'],data['stary'], 'o',color=data['dmag'])
    if data['dmag'] != None:
        if n.max(data['dmag']) != 0:
            pyl.figure()
            plot_contour_irregulargrid(data['starx'],data['stary'],data['dmag']*1000,
                                       nxbins=100, nybins=100, cb_title='mmag')
            pyl.title('Offsets applied across FoV (mmag)')
            pyl.savefig(fname,format='png')
    ui.sqlEndConnect(conn, cursor, verbose=True)
    

def nstarobs(table, fname='StatsNobs.png', file=None):
    """plot the distribution of the number of observations per star"""
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                               dbtype='pgsql', username=None, passwdname=None)
    cmd='select max(starid) from %s'%table
    cursor.execute(cmd)
    nid=cursor.fetchall()
    bins=n.arange(n.round(nid))+1
    cmd='select starid from %s'%table
    cursor.execute(cmd)
    temp=cursor.fetchmany(size = int(1e6))
    hist,bins=n.histogram(temp,bins=bins)
    bighist=hist+0
    temp=cursor.fetchmany(size = int(1e6))
    while (temp !=[]):
        hist,bins=n.histogram(temp,bins=bins)
        bighist=bighist+hist
        temp=cursor.fetchmany(size = int(1e6))
    #now I have an array that has how many observations each star had
    temp_num,temp_bins = n.histogram(bighist)
    med=n.median(bighist)
    #final_bins=n.arange(n.round(n.max(temp_bins)))+.5
    final_bins=n.arange(3.*med)+0.5
    pyl.figure()
    num,bins,patches = pyl.hist(bighist, bins=final_bins )
    pyl.title('Median=%.1f'%med)
    pyl.xlabel('Number of Repeat Observations')
    pyl.ylabel('Number of Stars')
    pyl.axvline(x=med)
    pyl.savefig(fname,format='png')
    if file != None:        
        n.savez(file, num=num,bins=bins, med=med)
    ui.sqlEndConnect(conn, cursor, verbose=True)

#taken from http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
def mapcount(filename):
    """read the number of lines in a file """
    f = open(filename, "r+")
    buf = mmap.mmap(f.fileno(), 0)
    lines = 0
    readline = buf.readline
    while readline():
        lines += 1
    return lines

def bigLoad(file, nentries=None, keys=None, dtype=float,skipline=0 ):
    """read a large file into a dictionary of numpy arrays """
    if nentries == None:
        nentries = mapcount(file)
        nentries = nentries - skipline
    #create large dictionary to hold data
    result={}
    for key in n.arange(n.size(keys)):
        result[keys[key]]=n.zeros(nentries, dtype=dtype[key])
    f = open(file, 'r')
    if skipline != 0:
        for i in range(skipline):
            f.readline()
    for i in range(nentries):
        line = f.readline()
        linevalues = line.split()
        for j in range(n.size(linevalues)):
            result[keys[j]][i]=linevalues[j]
    f.close()
    return result


def stupid_fitter():
    """Skip the fitter and just assume the zeropoints are known """    
    stardata = bigLoad('star_obs.dat.s', skipline=1, keys=('patchid', 'starid', 'starobsmag','magerr'), dtype=('int','int','float','float'))
    print('star_obs.dat loaded')
    patch_data = {}
    star_out = {}
    patch_data['upatch']=n.unique(stardata['patchid'])
    patch_data['patch_zp'] = n.zeros(n.size(patch_data['upatch']))
    star_out['ustar'] = n.unique(stardata['starid'])
    star_out['star_mags']=n.zeros(n.size(star_out['ustar']))
    print('setup the arrays, now to go through the sorted list')
    k=0
    for i in n.arange(n.size(star_out['ustar'])):
        good=n.array(k)
        while stardata['starid'][k+1] == stardata['starid'][k]:
            n.append(good,k)
            k += 1
        star_out['star_mags'][i]=n.average(stardata['starobsmag'][good],weights=stardata['magerr'][good])   
#    for i in n.arange(n.size(star_out['ustar'])):
#        good = n.where(stardata['starid'] == star_out['ustar'][i])
#        star_out['star_mags'][i]=n.average(stardata['starobsmag'][good],weights=stardata['magerr'][good])
    ui.writeDatafile('test_bestfit_Patch.dat', patch_data, ('upatch','patch_zp'), printheader=False)
    ui.writeDatafile('test_bestfit_Star.dat',star_out, ('ustar','star_mags'), printheader=False)
    
def rmsbins(name, file=None):
    """Find the RMS for different cloud levels"""
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                              dbtype='pgsql', username=None, passwdname=None)
    ranges=[0.,.3,.7,1.5]
    k=0
    cmd = 'SELECT AVG(f.starmagcalib- sd.startruemag) from '+name+'_bestfit as f, '+name+'_stars as sd where f.starid=sd.starid;'
    print(cmd)
    cursor.execute(cmd)
    mp=n.median(n.array(cursor.fetchall()))
    print('found floating zp as %f'%mp)
    for i in n.arange(n.size(ranges)-1):
        print(k)
        cmd = 'SELECT sd.starID, AVG(o.starObsMag - p.dmagcal - sd.startruemag), STDDEV(o.starObsMag - p.dmagcal - sd.startruemag) from '+name+'_bestfitpatch as p, '+\
        name+' as o, '+name+'_stars as sd, '+name+'_bestfit as sb WHERE p.patchid = o.PatchID and o.starid = sd.starID and o.starid = sb.starid and %f+p.dmagcal >= %f and %f+p.dmagcal < %f GROUP BY sd.starid;'%(mp,ranges[k], mp,ranges[k+1])
        print(cmd)
        cursor.execute(cmd)
        data = n.array(cursor.fetchall())
        if n.size(data) > 1:
            print('size of returned data = ', n.size(data), n.shape(data))
            floating_zp = n.median(data[:,1])
            data[:,1]=data[:,1]-floating_zp
            nstars=len(data[:,1])
            print('nstars =%i'%nstars)
            good=n.where(n.abs(data[:,1]) < 0.5)
            avs=n.ravel(data[good,1])
            stds=n.ravel(data[:,2])
            stds=stds[n.where((stds > 0) & (stds < .5))]
            del data
            pyl.figure()
            num,bins,patches = pyl.hist(avs*1000,bins=50)
            pyl.xlabel('<m$_{calibrated}$-m$_{true}$>$_i$ (mmag)')
            pyl.title('RMS=%.1f mmag, clipped %i stars'%(n.std(avs*1000),nstars-len(avs) ))
            pyl.figtext(.65,.75, '%.1f < Zp < %.1f'%(ranges[k],ranges[k+1]))
            pyl.savefig('Statsavebinstars%i.png'%(i),format='png')
        
            pyl.figure()
            num,bins,patches = pyl.hist(stds*1000,bins=50)
            if n.max(bins) > 150:
                pyl.figure()
                num,bins,patches = pyl.hist(stds*1000,bins=500, log=True)       
            pyl.xlabel('RMS(m$_{calibrated}$-m$_{true}$)$_i$ (mmag)')
            pyl.title('Repeatability: Median=%.1f mmag, clipped %i stars' %(n.median(stds*1000), nstars-len(stds)))#, ranges[k], ranges[k+1])
            pyl.figtext(.65,.75, '%.1f < Zp < %.1f'%(ranges[k],ranges[k+1]))
            pyl.savefig('StatsRMSSbins%i.png'%i,format='png')
        k=k+1
    ui.sqlEndConnect(conn, cursor, verbose=True)

def rmsrad(name, bins=10):
    """plot the RMS of the residuals as a function of radius on the FoV"""
    conn,cursor=ui.sqlConnect('localhost', dbname='calsim',\
                              dbtype='pgsql', username=None, passwdname=None)

    cmd = 'SELECT AVG(f.starmagcalib- sd.startruemag) from '+name+'_bestfit as f, '+name+'_stars as sd where f.starid=sd.starid;'
    print(cmd)
    cursor.execute(cmd)
    mp=n.median(n.array(cursor.fetchall()))
    print('found floating zp as %f'%mp)

    cmd="SELECT max( SQRT(starx^2+stary^2)) from "+name+";"
    cursor.execute(cmd)
    rmax=n.ravel(n.array(cursor.fetchall()))
    radii=n.arange(bins+1.)/bins*rmax
    mean_rms=n.zeros(bins)
    spread_rms=n.zeros(bins)
    for i in n.arange(bins):
        cmd = "SELECT AVG(o.starObsMag - p.dmagcal - sd.startruemag) from "+name+"_bestfitpatch as p, "+\
        name+" as o, "+name+"_stars as sd, "+name+"_bestfit as sb WHERE p.patchid = o.PatchID and"+\
        " o.starid = sd.starID and o.starid = sb.starid and sqrt(o.starx^2+o.stary^2) between %f AND %f GROUP BY sd.starid;"%(radii[i],radii[i+1])
        cursor.execute(cmd)
        data=n.ravel(n.array(cursor.fetchall()))
        #print 'data =',data
        good=n.where(n.abs(data-mp) < 0.1) #clip any massive outliers
        mean_rms[i]=n.std(data[good])
        #spread_rms[i]=n.std(data)
        #print 'min =%f, max=%f'%(n.min(data), n.max(data))
    pyl.figure()
    #print mean_rms
    #print spread_rms
    #pyl.errorbar(radii[1:bins+1]/rmax, mean_rms*1000, yerr=spread_rms*1000, fmt='o', )
    pyl.plot(radii[1:bins+1]/rmax, mean_rms*1000, 'ko')
    pyl.xlim([0,1.1])
    pyl.xlabel('Radius (fraction of FoV)')
    pyl.ylabel('Residual RMS (mmag)')
    pyl.savefig('StatsRMSrad.png',format='png')
    ui.sqlEndConnect(conn, cursor, verbose=True)

