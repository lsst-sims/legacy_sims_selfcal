from __future__ import print_function
#####
#  Lynne Jones, ljones@astro.washington.edu.
#
# This class has some input / output methods for python which are useful for
# general analysis. Suitable for flexible, quick input/output, but if you need a
# large data volume or faster read/write, you might want something more efficient.
#
# This basically, gives you a set of functions to connect to, then query, then
# close a database or read from a data file,
# and place the data read from those sources into a python dictionary of
# numpy arrays (keys for the dictionary are user-defined).
# 
# Requires sqlalchemy module for python; depending on which type of database you connect to, may
#  then require MySQLdb or pgdb, etc. (sqlite is built into python).
#
# * sqlConnect(hostname='localhost', username='lsst', passwdname='lsst', dbname='opsimdev',
#             type="mysql", verbose=True):
#   This sets up your connection to the database. Returns the connection and cursor objects.
# * sqlEndConnect(conn, cursor, verbose=True):
#   This ends your connection to the database. Execute when done with all queries.
# * sqlQuery(cursor, query, verbose=True):
#   This executes your query.
# * sqlQueryLimited(cursor, query, verbose=True, offset=0, drows=1000):
#   This executes your query, for a subset of the data (drows). Might not work for postgres?
# * assignResults(sqlresults, keys):
#   This puts your database results into the python dictionary of numpy arrays.
# * readDatafile(infilename, keys, keytypes=None,
#   This reads your data file and puts the results into numpy arrays. If the
#    file includes non-numeric data, include the keytypes.
# * writeDatafile(outfilename, data, keys, printheader=True, newfile=True):
#   Take your python dictionary and write it out to disk.
#
#####

# general python modules
import numpy
import sqlalchemy
from sqlalchemy import create_engine

# for opsim
minexpmjd = 49353.
maxexpmjd = 53003.
midnight = 0.666


class InputOutput:

    def __init__(self):
        """Initialize input/output class. 
            This class simply keeps track of sql connection values and provides similar methods for data files."""
        pass 
        return

    ## Methods to deal with sql connection, query, sending results.
    
    def sqlConnect(self, hostname='localhost', username='lsst', password='lsst', database='opsimdev', 
                   dbtype="mysql", driver=None, verbose=False):
        """Connect to database."""
        # The sqlalchemy URI is dialect(dbtype)+driver://username:password@host:port/database
        # Examples: dialect/drivers 'mssql+pyodbc://...' , 'sqlite' (default driver=pysqlite)
        #    and 'mysql' (default driver=MySQLdb) , 'postgresql' (default driver=psycopg2)
        self.hostname = hostname
        self.username = username
        self.password = password
        self.database = database
        self.dbtype = dbtype
        self.driver = driver        
        if ((self.dbtype=='sqlserver') & (self.username=='lsst')):
            # Use default LSST username/password to mssql server.
            self.hostname = 'SQLSERVERDB'
            self.username = 'LSST-2'
            self.password = 'L$$TUser'
            self.dbtype = 'mssql'
            self.driver = 'pyodbc'
            self.engine = create_engine('mssql+pyodbc://LSST-2:L$$TUser@SQLSERVERDB', echo=verbose)
        # sqlite has some slightly different formatting. 
        elif (self.dbtype=='sqlite'):
            if self.database!='memory':
                # Assume that this is a real file that should be opened as an sqlite table,
                #  either for reading or writing.
                self.engine = create_engine('sqlite:///%s' %(self.database))
            else:
                # Create a new in-memory table. 
                self.engine = create_engine('sqlite://')
        else:
            # Any other kind of connection -  mysql or postgres, most likely.
            if self.driver==None:
                self.engine = create_engine('%s://%s:%s@%s/%s' %(self.dbtype, self.username, self.password,
                                                                  self.hostname, self.database), echo=verbose)    
            else:
                self.engine = create_engine('%s+%s://%s:%s@%s/%s' %(self.dbtype, self.driver, self.username, self.password, 
                                                                     self.hostname, self.database), echo=verbose)   
        if verbose:
            print("Connected to %s database %s as %s." %(self.dbtype, self.database, self.username))    
        self.conn = self.engine.connect()
        return 

    def sqlEndConnect(self, verbose=True):
        """End connection to database."""
        try:
            self.conn.close()
        except AttributeError:
            if verbose:
                print('# No connection existed to close.')
        if verbose:
            print('# Closed connection to database.')
        return

    def sqlQuery(self, query, verbose=True):
        """Send query to db table, return results."""
        if verbose:
            print('#', query)
        # Execute query.
        resultproxy = self.conn.execute(query)
        sqlresults = resultproxy.fetchall()
        resultproxy.close()
        return sqlresults

    def sqlCreateTable(self, table, keys, keytypes, stringlength=10, verbose=True):
        """Create table for the db specified by self.conn. 
            keys = column names, keytypes = column types."""
        # Check for 'string' or 'str' in keytypes: this has to be changed to CHAR or VARCHAR.
        for kt in keytypes:
            if ((kt=='string') | (kt=='str')):
                kt = 'char(%d)' %(stringlength)
                if verbose:
                    print('#! Warning - heads up, using strings of length %d' %(stringlength))
        query = 'create table %s (' %(table)    
        for i in range(len(keys)):  
            query = query + '%s %s, ' %(keys[i], keytypes[i])
        # Remove the last comma.
        query = query[:-2]
        # Close the parenthesis.
        query = query + ')'
        if verbose:
            print('# %s' %(query))
        # Begin transaction. 
        trans = self.conn.begin()
        self.conn.execute(query)
        trans.commit()
        return 

    def sqlWriteData(self, table, data, keys, data_rows=False, verbose=True):
        """Write data (dict. of numpy arrays) into a database specified by self.conn, where keys = columns. 
        If data_rows is true, data is expected to be [(a1 b1 c1), (a2, b2, c2)] format. """
        basequery = 'insert into %s (' %(table)
        for k in keys:
            basequery = basequery + '%s, ' %(k)
        basequery = basequery[:-2]
        basequery = basequery + ') values '
        trans = self.conn.begin()
        # If data given as 'sqlresults' style (
        if data_rows:
            data_len = len(data)
            for row in data:
                query = basequery + '%s' %(row)
                self.conn.execute(query)
                # At end of updating of table, commit / flush to table.
        # Else data given as dictionary of numpy arrays.
        else:
            basequery = basequery + '('
            data_len = len(data[keys[0]])
            for i in range(data_len):
                for k in keys:
                    query = basequery + '%s, ' %(data[k][i])
                query = query[:-2]
                query = query + ')'
                self.conn.execute(query)
        trans.commit()
        if verbose:
            print("# Updated table %s with %s rows of data" %(table, data_len))
        return

    def sqlCreateIndexes(self, table, keys_to_index, verbose=True):
        """Create indexes on keys_to_index in table. """
        for k in keys_to_index:
            query = 'create index %s_idx on %s (%s)' %(k, table, k)
            if verbose:
                print('# %s' %(query))
            self.conn.execute(query)
        return

    def assignResults(self, sqlresults, keys, stringlength=10, verbose=True):
        """Split sqlresults based on keys, return a dictionary of numpy arrays with results.
            Keys must map to the sqlquery results. """
        # Note that this function requires making a copy of all of your data, so 
        # may be less than ideal for very large data sets. 
        arraylen = len(sqlresults)
        # Return empty dictionary if number of results were 0.
        if arraylen==0:
            value = {}
            for key in keys:
                value[key] = None
            return value
        # Test data types of the sqlresults.
        result = sqlresults[0]
        j = 0
        datatype = {}
        for key in keys:
            test = result[j] 
            j = j+1
            datatype[key] = type(test)
        # Set up the dictionary for the data values.
        value = {}
        for key in keys:
            if datatype[key] == str:
                value[key] = numpy.empty((arraylen,), dtype=(datatype[key], stringlength))
                if verbose:
                    print('#! Warning - heads up - using strings of length %d' %(stringlength))
            else:
                value[key] = numpy.empty((arraylen,), dtype=datatype[key])
        # Start assigning sql query results.
        i = 0
        for result in sqlresults: 
            j = 0  
            for key in keys:
                value[key][i] = result[j]
                j = j+1
            i = i+1
        # Done assigning sql query results.
        return value

    ## Crossover: read a sqlite database on disk into an in-memory database.

    def sqlDiskToMem_DB(self, filename_database, tables):
        """Read on-disk sqlite database file into memory, copying all tables specified in 'tables'. """
        if type(tables) == str:
            tables = [tables,]
        # End existing connection. 
        self.sqlEndConnect()
        # Set up a connection to an inmemory sqlite db.        
        self.sqlConnect(dbtype='sqlite', database='memory')
        # Note that database_filename must be the FULL path to the file. 
        self.conn.execute('attach "%s" as diskdb' %(filename_database))
        for t in tables:
            self.conn.execute('create table %s as select * from diskdb.%s' %(t, t))
        self.conn.execute('detach diskdb')
        return

    def sqlDiskToMem(self, filename_flat, table, keys, keytypes=None):
        """Read on-disk flat file into memory, creating a table with columns keys of keytypes."""
        raise Exception("Not implemented yet")
        return
        
        
    ## Routines to deal with on-disk files with similar methods.

    def readDatafile(self, infilename, keys, keytypes=None,
                     startcol=0, sample=1, skiprows=0, delimiter=None, verbose=True):
        """Read values from infilename, in columns keys, and return dictionary of numpy arrays.
    
            Limitation - while startcol can be >0, cannot skip columns beyond that point."""
        # Attempt to open file.
        f = open(infilename, 'r')
        if verbose:
            print("# Reading file %s" %(infilename))
        # Read data from file.
        value = {}
        for key in keys:
            value[key] = []
        sampleskip = 0
        for line in f:
            # Ignore line if it's obviously a comment.
            if (line.startswith('!') | line.startswith('#')):
                continue
            # Increment counter on whether to read/assign the values from this line or not.
            sampleskip = sampleskip + 1
            # Assign values for each line we want to keep. 
            if sampleskip == sample:
                linevalues = line.split()
                if delimiter!=None:
                    linevalues = line.split()
                # j is a column counter for linevalues.
                j = startcol
                # If the line is long enough, assign values (this will skip incomplete lines without errors). 
                if (len(linevalues)-j)>=len(keys):
                    # Loop through columns, assign to keys. 
                    for key in keys:
                        value[key].append(linevalues[j])
                        j = j+1
                sampleskip = 0
        # End of loop over lines.
        # Convert each list to numpy arrays.
        for key in keys:            
            if keytypes!=None:
                value[key] = numpy.array(value[key], dtype=keytypes[keys.index(key)])
            else:
                value[key] = numpy.array(value[key], dtype='float')
        f.close()
        return value

    def writeDatafile(self, outfilename, data, keys, printheader=True, newfile=True):
        """Write dictionary of numpy arrays to file, potentially appending to existing file."""
        import sys
        if newfile:
            try:
                fout = open(outfilename, 'w')
            except IOError:
                print("Couldn't open new file %s" %(outfilename), file=sys.stderr)
                return
        else:
            try:
                fout = open(outfilename, 'a')
            except IOError:
                print("Couldn't open file %s for appending." %(outfilename), file=sys.stderr)
                return
        # Print header information if desired.
        if printheader: 
            writestring = "# "
            for key in keys:
                writestring = writestring + " %s" %(key)
            print(writestring, file=fout)
        # Print data information.
        for i in range(len(data[keys[0]])):
            writestring = ""
            for key in keys:
                writestring = writestring + "%s " %(data[key][i])
            # Trim trailing white space.
            writestring = writestring[:-1]
            # Done constructing write line.
            print(writestring, file=fout)
        # Done writing data.
        fout.close()
        return 
