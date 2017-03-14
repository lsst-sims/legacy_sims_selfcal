from __future__ import print_function
#Take the output files from simSelfCalib.py and send them over to a Postgres database
#assumes you are running the script in the directory where files are to be coppied from
#
#
#
#
#
#
#
#
#


import numpy as n
import pgdb as pg
import os

def fitted2db(master_table_name, index_tag='_', dbname='localhost:calsim', bestfitfile='test_bestfit_Star.dat', bestfitpatchfile='test_bestfit_Patch.dat'):
    """Copy the results from the selfcal fitter to a database """
    bestfit_table_name = master_table_name+'_bestfit'
    bestfit_patch_table_name = master_table_name+'_bestfitpatch'

    #reformat simple.x output    
    os.system("perl -i -p -e 's/^\s+//;' "+bestfitfile) ##remove leading spaces
    os.system("perl -i -p -e 'tr/ //s;' "+bestfitfile)  ##remove  any extra spaces
    os.system("perl -i -p -e 's/^\s+//;' "+bestfitpatchfile) ##remove leading spaces
    os.system("perl -i -p -e 'tr/ //s;' "+bestfitpatchfile)  ##remove  any extra spaces

    #open the database
    print('opening connection')
    db = pg.connect(dsn = dbname)
    cursor = db.cursor ()
    cwd=os.getcwd()

    #Create the bestfit stars table
    cursor.execute ("DROP TABLE IF EXISTS "+bestfit_table_name)
    print("DROP TABLE IF EXISTS "+bestfit_table_name)
    cmd = "CREATE TABLE " + bestfit_table_name + " (StarID int, StarMagCalib real);"
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
    cmd = "CREATE TABLE " + bestfit_patch_table_name + " (PatchID int, dmagcal real);"
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
    cursor.execute('DROP INDEX IF EXISTS master_patchid'+index_tag+' ;')

    cmd = "CREATE INDEX patchid_patch" +index_tag +" on " + bestfit_patch_table_name + " (PatchID);"
    cursor.execute (cmd)

    cmd = "CREATE INDEX starid_bf" + index_tag +" on " +  bestfit_table_name + " (StarID);"
    cursor.execute (cmd)
    
    db.commit()
    print('finished, closing connection')
    cursor.close()
    db.close()



def sim2db(master_table_name, index_tag='_', dbname='localhost:calsim', masterfile='master_cal.dat', starfile='stardata.dat',
           bestfitfile='test_bestfit_Star.dat', bestfitpatchfile='test_bestfit_Patch.dat', 
           star_obs=None):
    """Take the simulation output from master_table.dat and the bestfit files. 
    Put them into a postgres database.  Index the tables on useful quanities (like StarID).  """

    #open the database
    print('opening connection')
    db = pg.connect(dsn = dbname)
    cursor = db.cursor ()
    cwd=os.getcwd()

    #copy in the stardata
    cursor.execute("DROP TABLE IF EXISTS "+master_table_name+"_stars")
    cmd = "CREATE TABLE "+master_table_name+ \
          "_stars (StarID int, StarDBID int, RA real, Dec real," + \
          " StarTrueMag real, StarColorTrue real);"
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
        StarObsMag  real, \
        StarMagErr  real, \
        dMag_var	 real, \
        dMag_snr	 real, \
        dMag_zp	 real, \
        dMag_zp_dist real, \
        dMag_color	  real, \
        dMag_color_dist  real, \
        dMag_sin	      real, \
        dMag_flat_sin	      real, \
        dMag_cloud            real, \
        dMag_cloud_image     real, \
        dMag_rr           real, \
        StarTrueMag      real, \
        StarTrueCol      real, \
        StarRA	      double precision, \
        StarDec	      double precision, \
        RAcen_fov	      real, \
        DecCen_fov	      real, \
        StarX	      real, \
        StarY	      real,\
        night         real,\
        dm5           real,\
        dMag_kep      real\
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

        
    db.commit()
    print('finished, closing connection')
    cursor.close()
    db.close()

########################
if __name__ == "__main__":
    # Main executable. 

    import sys
    if len(sys.argv) < 2:
        dbname = 'localhost:calsim'
    else:
        dbname = sys.argv[1]
    if len(sys.argv) < 3:
        simname = 'selfcal_v01'
    else:
        simname = sys.argv[2]
    sim2db(simname, dbname = dbname)
    fitted2db(simname, dbname = dbname)
    
