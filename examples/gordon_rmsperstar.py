import numpy as np
import time 
import pylab as pyl

#attempt # 2 at optimizing a large file read


def readDatafile(infilename, keys, dtype=None,
                 startcol=0, cols=None, sample=1, skiprows=0, delimiter=None):
    """Read values from infilename, in columns keys, and return structured numpy array.
    
    Limitation - while startcol can be >0, cannot skip columns beyond that point."""
    # this can be done with loadtxt too, if you don't unpack (put in recarray) -- TODO
    # open file
    import sys
    try:
        f = open(infilename, 'r')
    except IOError:
        print >>sys.stderr, "Couldn't open file %s" %(infilename)
        print >>sys.stderr, "Returning None values"
        value = {}
        for i in keys:
            value[i] = None
        return value
    print >>sys.stderr, "Reading file %s" %(infilename)
    # Read data from file.
    value = {}
    for key in keys:
        value[key] = []
    sampleskip = 0
    i = 0
    for line in f:
        i = i + 1
        if line.startswith("!"):
            continue
        if line.startswith("#"):
            continue
        sampleskip = sampleskip + 1
        # Assign values
        if sampleskip == sample:
            linevalues = line.split()
            j = startcol
            if len(linevalues)>=(len(keys)-j):
                for key in keys:
                    try:
                        value[key].append(linevalues[j])
                    except IndexError:
                        print "readDataFile failed at line %d, column %d, key=%s" \
                              %(i, j+1, key)
                        print "Data values: %s" %(linevalues)
                        raise IndexError
                    j = j+1
            sampleskip = 0
    # End of loop over lines.
    # Convert dictionary to numpy array
    result = np.zeros(np.size(value[keys[0]]),dtype=dtype)
    for name in keys:
        result[name]=np.array(value[name])
        del value[name]
    #for key in keys:
    #    if keytypes!=None:
    #        value[key] = np.array(value[key], dtype=keytypes[keys.index(key)])
    #    else:
    #        value[key] = np.array(value[key], dtype='float')
    f.close()
    return result


t1=time.time()
patchfit=readDatafile('test_bestfit_Patch.dat', ['patchid','zp'] ,
                     dtype=[('patchid', '<i8'),  ('zp', '<f8')], skiprows=1)
t2=time.time()
print 'time to read patch', t2-t1
t11=t1

t1=time.time()
starobs=readDatafile('star_obs.dat', ['patchid','starid','obsmag','magerr'] ,
                     dtype=[('patchid', '<i8'), ('starid', '<i8'), ('obsmag', '<f8'), ('magerr', '<f8')], skiprows=1)
t2=time.time()
print 'time to read starobs.dat', t2-t1

t1=time.time()
starobs=np.sort(starobs,order='patchid') #order by patch id
#can try changing the sort type to look for speedup--nope, need to search for more, maybe parallel or pre-sort
t2=time.time()
print 'time to sort by patch id=', t2-t1

print 'starting searchsorted left'
pl=np.searchsorted(starobs['patchid'], patchfit['patchid'])
print 'starting search sorted right'
pr=np.searchsorted(starobs['patchid'], patchfit['patchid'], side='right')
#might need to do a trim patches from starobs

print 'constructing zeropoint array'
zparray=starobs['obsmag']*0
for i in xrange(np.size(pl)):
    if pl[i] != pr[i]:
        zparray[pl[i]:pr[i]]=zparray[pl[i]:pr[i]]+patchfit['zp'][i]

print 'subtracting zeropoint array'
starobs['obsmag']=starobs['obsmag']-zparray

print 'sorting by starid'
t1=time.time()
starobs.sort(order='starid')
t2=time.time()
print 'time to sort by starid =',t2-t1


ids=np.unique(starobs['starid'])
stds=np.zeros(np.size(ids))
print 'starting searchsorted left'
left=np.searchsorted(starobs['starid'], ids)
print 'starting search sorted right'
right=np.searchsorted(starobs['starid'], ids,'right')

print 'computing stdev for each star'
for i in xrange(np.size(ids)):
    stds[i]=starobs['obsmag'][left[i]:right[i]].std()

print 'median RMS per star', np.median(stds)

pyl.figure()
num,bins,patches = pyl.hist(stds*1000,bins=50)
if np.max(bins) > 150:
    pyl.figure()
    num,bins,patches = pyl.hist(stds*1000,bins=500, log=True)       
pyl.xlabel('RMS(m$_{calibrated}$)$_i$ (mmag)')
pyl.title('Repeatability: Median=%.1f mmag'
              %(np.median(stds*1000)))
pyl.savefig('rmsperstar.png', format='png')

t2=time.time()
print 'total runtime = %.2f min'%((t11-t2)/60)


#now to import pylab and make the histogram

