#####
#  Lynne Jones, ljones@astro.washington.edu
#  svn version info : $Id$
#
#
#####

# general python modules
import numpy as n
import pylab as pyl
# self-calibration modules
import lsst.sims.selfcal.analysis.useful_input as ui

#calibtable = "bhb"
#calibtable = "msrgb"

deg2rad = n.pi/180.0
rad2deg = 180.0/n.pi


def connectDB():
    """Connect to the default calsim database, using parameters from environment variables.
    Returns connection and cursor. """
    calibhost, calibuser, calibpasswd, calibdb, calibdbtype = ui.setupConnEnv("LSST_CALSIM_")
    conn, cursor = ui.sqlConnect(hostname=calibhost, username=calibuser,
                                 passwdname=calibpasswd, dbname=calibdb,
                                 type=calibdbtype, verbose=True)
    return conn, cursor

def endConnectDB(conn, cursor):
    ui.sqlEndConnect(conn, cursor)
    return

def getData(cursor, calibtable, queryrestriction=None,
            querykeys=('ra', 'decl', 'rmag'), labelkeys=None):
    """Sample or basic query to the calsim database. Returns a dictionary of numpy arrays."""
    if (type(querykeys).__name__ != 'tuple') & (type(querykeys).__name__ != 'list'):
        raise Exception("Expecting querykeys to be a tuple or list.")
    if labelkeys == None:
        labelkeys = querykeys
    if len(labelkeys)!=len(querykeys):
        raise Exception("Querykeys and labelkeys must be the same length.")
    # Get on with things.
    sqlquery = "select %s" %(querykeys[0])
    for i in range(1, len(querykeys)):
        sqlquery = sqlquery + ", %s" %(querykeys[i])
    sqlquery = sqlquery + " from %s" %(calibtable)
    if queryrestriction!=None:
        sqlquery = sqlquery + " %s" %(queryrestriction)
    sqlresults = ui.sqlQuery(cursor, sqlquery)
    data = ui.assignResults(sqlresults, querykeys)
    return data

def wrapRA(ra):
    """Wrap RA into range of 0 to 360 degrees. Uses radians."""
    ra = ra % (360.0*deg2rad)
    return ra

def capDec(dec):
    """Cap declination to -90 to 90 degrees. Does not wrap, just truncates. Uses radians."""
    dec = dec.clip(-90*deg2rad, 90*deg2rad)
    return dec

def wrapLon(lon, lonCen):
    """Wrap longitude into -180 to 180 instead of 0 to 360. Uses radians. """
    lontemp = lon - lonCen
    lontemp = lontemp % (360.0*deg2rad)
    condition = (lontemp > (180*deg2rad))
    lontemp[condition ] = lontemp[condition] - (360.0*deg2rad)
    return lontemp

def hammer_project_toxy(ra, dec, raCen=0, decCen=0):
    """Calculate x/y projection of lon/lat (or RA/Dec) in a hammer projection.
    Nature of projection means that decCen always = 0.
    Input radians. Returns x/y. """
    if (ra.shape[0]!=dec.shape[0]):
        raise Exception("Expect lon and lat input to hammer_project_toxy to be same length.")
    if decCen!=0:
        print "Ignoring decCen=0; this is a hammer projection (decCen=0 automatically)."
    rt2=n.sqrt(2.0)
    # wrap longitude so lonCen is at center, and appropriately phased for cosine
    # (i.e. lon 0-360 must go to -180 to 180, with 0 at center.)
    ratemp = wrapLon(ra, raCen)
    denom=n.sqrt(1.0+(n.cos(dec)*n.cos(ratemp/2.0)))
    x = (2.0*rt2*n.cos(dec)*n.sin(ratemp/2.0))/denom
    y = rt2*n.sin(dec)/denom
    return x,y

def hammer_project_tosky(x, y, raCen=0, decCen=0):
    """Calculate sky ra/dec from hammer projection x/y coordinates.
    Nature of the projection means that decCen=0 always. Input raCen in radians. 
    Returns ra/dec. """
    if (x.shape[0]!=y.shape[0]):
        raise Exception("Expect x and y to hammer_project_tosky to be the same length.")
    if decCen!=0:
        print "Ignoring decCen, setting =0 as this is a hammer projection."
    z = n.sqrt(1 - (x/4.0)**2 - (y/2.0)**2)
    ra = 2 * n.arctan2(z*x, (2*(2*z**2 - 1))) + raCen
    ra = wrapRA(ra)
    dec = n.arcsin(z*y)
    return ra, dec

def gnomonic_project_toxy(ra, dec, raCen, decCen):
    """Calculate x/y projection of RA1/Dec1 in system with center at RAcen, Deccen.
    Input radians. Returns x/y."""
    # also used in Global Telescope Network website
    if (len(ra)!=len(dec)):
        raise Exception("Expect RA and Dec arrays input to gnomonic projection to be same length.")
    cosc = n.sin(decCen) * n.sin(dec) + n.cos(decCen) * n.cos(dec) * n.cos(ra-raCen)
    x = n.cos(dec) * n.sin(ra-raCen) / cosc
    y = (n.cos(decCen)*n.sin(dec) - n.sin(decCen)*n.cos(dec)*n.cos(ra-raCen)) / cosc
    return x, y

def gnomonic_project_tosky(x, y, raCen, decCen):
    """Calculate RA/Dec on sky of object with x/y and RA/Cen of field of view.
    Returns Ra/Dec in radians."""
    denom = n.cos(decCen) - y * n.sin(decCen)
    ra = raCen + n.arctan2(x, denom)
    dec = n.arctan2(n.sin(decCen) + y * n.cos(decCen), n.sqrt(x*x + denom*denom))
    return ra, dec

def quick_hammer(stardat, raCen=0):
    stardat['x'], stardat['y'] = hammer_project_toxy(stardat['ra']*deg2rad, stardat['dec']*deg2rad,
                                                     raCen*deg2rad)
    return stardat

def make_radec_grid(projectionmethod, raCen=0, decCen=0, rastep=20, decstep=20,
                    ralabel_dec=0, declabel_ra=0, ralabel_step=2, declabel_step=1,
                    xgridlim=None, ygridlim=None, newfig=False):
    """Add an ra/dec grid, onto a particular projection.
    Give this the center pointing for Ra(optional: Dec) in radians, optional unless !=0. 
    The rastep / decstep for grid lines are in degrees, optional.
    The ralabel_dec / declabel_ra for ra/dec label locations are in degrees, optional.
    x/ygridlim are the x/y limits (in xy units) of the final plot, optional. """
    if newfig:
        pyl.figure()
    pyl.axis('equal')
    # Add ra lines. 
    dec = n.arange(-90, 91, 1, dtype='float')
    dec = dec * deg2rad
    for r in range(0, 360, rastep):
        ra = n.zeros(len(dec), dtype='float')
        ra = ra + r
        ra = ra * deg2rad
        x, y = projectionmethod(ra, dec, raCen, decCen)
        pyl.plot(x, y, 'k:')
    # Make sure have lines at both edges of plot.
    redge = [(raCen + 179.9*deg2rad), (raCen - 179.9*deg2rad)]
    for r in redge:
        ra = n.zeros(len(dec), dtype='float')
        ra = ra + r
        x, y = projectionmethod(ra, dec, raCen, decCen)
        pyl.plot(x, y, 'k-')
    # Label RA lines
    # Use xgridlim/ygridlim to remove text which will fall outside final boundaries of plot.
    raoffset = 1
    decoffset = decstep/4.0 
    ra = n.arange(0, 360, rastep*ralabel_step) + raoffset
    dec = n.zeros(len(ra), dtype='float') + decoffset + ralabel_dec
    x, y = projectionmethod(ra*deg2rad, dec*deg2rad, raCen, decCen)
    ra = ra - raoffset
    for i in range(len(ra)):
        plot_this_text = True
        if xgridlim!=None:
            if (x[i]>=xgridlim[1]) | (x[i]<=xgridlim[0]):
                plot_this_text = False
        if ygridlim!=None:
            if (y[i]>=ygridlim[1]) | (y[i]<=ygridlim[0]):
                plot_this_text = False
        if plot_this_text:
            pyl.text(x[i], y[i], "%.0f" %(ra[i]))
    # Add dec lines.
    ra = n.arange(0, 360, 1, dtype='float')
    ra = ra * deg2rad
    for d in range(-80, 90, decstep):
        dec = n.zeros(len(ra), dtype='float')
        dec = dec + d
        dec = dec * deg2rad
        x, y = projectionmethod(ra, dec, 0, decCen)
        # Split x for pretty plotting.
        condition = (x > 0)
        x1 = x[condition]
        y1 = y[condition]
        condition = (x < -0.001 )
        x2 = x[condition]
        y2 = y[condition]
        pyl.plot(x1, y1, 'k:')
        pyl.plot(x2, y2, 'k:')
    # Label dec lines
    raoffset = rastep/2
    decoffset = 1
    dec = n.arange(-80, 90, decstep*declabel_step) + decoffset
    ra = n.zeros(len(dec), dtype='float') + raoffset + declabel_ra
    x, y = projectionmethod(ra*deg2rad, dec*deg2rad, raCen, decCen)
    dec = dec - decoffset
    for i in range(len(ra)):
        if dec[i]!=0:
            plot_this_text = True
            if xgridlim!=None:
                if (x[i]>=xgridlim[1]) | (x[i]<=xgridlim[0]):
                    plot_this_text = False
            if ygridlim!=None:
                if (y[i]>=ygridlim[1]) | (y[i]<=ygridlim[0]):
                    plot_this_text = False
            if plot_this_text:
                pyl.text(x[i], y[i], "%.0f" %(dec[i]))
    # Done.
    return
            
def calc_gridsize(data=None, binsize='fov', rad_fov=1.8, patches_per_side=4):
    """Calculate number of bins for pyl.hexbin plots, appropriate for either the
    full field of view (binsize=fov) or a single patch (binsize=patch).
    Calculation uses data['ra'], data['dec'], in degrees, as well as radius_fov (1.8)
    and patches_per_side (4).
    Returns gridsize, =[gridsizex, gridsizey]."""
    binsize_opts = ['fov', 'patch']
    if binsize not in binsize_opts:
        raise Exception("binsize must be one of %s" %(binsize_opts))
    sqrt2 = n.sqrt(2)
    if binsize=='fov':
        bin_len = rad_fov*2/sqrt2
    elif binsize=='patch':
        bin_len = rad_fov*2/sqrt2/patches_per_side
    if data != None:
        gridsizex = round((data['ra'].max() - data['ra'].min())/(bin_len))
        gridsizey = round((data['dec'].max() - data['dec'].min())/(bin_len))
    else:
        gridsizex = round(360. / bin_len)
        gridsizey = round(90. / bin_len)
    gridsize = [int(gridsizex), int(gridsizey)]
    print 'gridsizes=', gridsize
    return gridsize


def calc_squaredeg(projection=hammer_project_toxy, size=1.0):
    """Calculate the x/y projection length of a patch with side of size 'size' (deg)."""
    ra = n.array([0.0, size*deg2rad])
    dec = n.array([0.0, 0.0])
    x, y = projection(ra, dec)
    dx = x[1]-x[0]
    ra = n.array([0.0, 0.0])
    dec = n.array([0.0, size*deg2rad])
    x, y = projection(ra, dec)
    dy = y[1]-y[0]
    return [dx, dy]


def plot_density(x, y, z=None, z_method=n.mean, zlimits=None, zscale=None,
                 gridsize=[360, 180], projectionmethod=hammer_project_toxy,
                 radecgrid=True, raCen=0,
                 xgridlim=None, ygridlim=None, newfig=True, cb_label=''):
    """Make a density plot of x/y, potentially with z of same length. Adds radec grid
    in hammer_projection. gridsize specifies binning size within plot.
    z_method (reduce_C_method), zlimits (min/max of z), and zscale (None/'log'/integer/sequence)
    affect plotting of z value, if provided.
    xgridlim, ygridlim specify limits on final plot (may be non-functioning). """
    cmap = None
    if newfig:
        pyl.figure()
    if z!=None:
        if zlimits == None:
            hx = pyl.hexbin(x, y, C=z, reduce_C_function=z_method, gridsize=gridsize,
                            bins=zscale, cmap=cmap)        
        else:
            delta = 0.00001
            hx = pyl.hexbin(x, y, C=z, reduce_C_function=z_method, gridsize=gridsize,
                            vmin=zlimits[0], vmax=zlimits[1]+delta, bins=zscale, cmap=cmap)
    else:
        if zlimits == None:
            hx = pyl.hexbin(x, y, gridsize=gridsize, mincnt=0.5, cmap=cmap, bins=zscale)
        else:
            hx = pyl.hexbin(x, y, gridsize=gridsize, mincnt=0.5,
                            vmin=zlimits[0], vmax=zlimits[1], cmap=cmap, bins=zscale)
    if radecgrid:
        make_radec_grid(projectionmethod, raCen=raCen*deg2rad,
                        xgridlim=xgridlim, ygridlim=ygridlim, newfig=False)
    pyl.xlabel("x")
    pyl.ylabel("y")
    if zlimits!=None:
        dtick = (zlimits[1] - zlimits[0])/5.0
        colorticks = n.arange(zlimits[0], zlimits[1]+dtick, dtick)
        cb = pyl.colorbar(hx, orientation='horizontal', pad=.1, fraction=0.05, shrink=1, aspect=40,
                          ticks=colorticks)
    else:
        cb = pyl.colorbar(hx, orientation='horizontal', pad=.1, fraction=0.05, shrink=1, aspect=40)
    if zscale == 'log':
        cb.set_label('log$_{10}$(N)')
    cb.set_label(cb_label)
    return 

def plot_contour_irregulargrid(x, y, z, nxbins=360, nybins=180):
    """Make a contour plot of irregularly gridded data."""
    # Grid the data.
    dx = (max(x) - min(x))/float(nxbins)
    dy = (max(y) - min(y))/float(nybins)
    xi = n.arange(min(x), max(x), dx)
    yi = n.arange(min(y), max(y), dy)
    xi, yi = n.meshgrid(xi, yi)
    zi = pyl.griddata(x, y, z, xi, yi)
    # Contour the gridded data.
    CS = pyl.contour(xi, yi, zi)
    CS = pyl.contourf(xi, yi, zi)
    pyl.colorbar()
    # Plot data points
    #pyl.scatter(x, y, 'b.')
    return

def count_number(x, y, xbinsize=None, ybinsize=None, nxbins=None, nybins=None):
    # Set up grid for contour/density plot.
    xmin = min(x)
    ymin = min(y)
    if (xbinsize!=None) & (ybinsize!=None):
        xbins = n.arange(xmin, max(x), xbinsize)
        ybins = n.arange(ymin, max(y), ybinsize)
        nxbins = xbins.shape[0]
        nybins = ybins.shape[0]
    elif (nxbins!=None) & (nybins!=None):
        xbinsize = (max(x) - xmin)/float(nxbins)
        ybinsize = (max(y) - ymin)/float(nybins)
        xbins = n.arange(xmin, max(x), xbinsize)
        ybins = n.arange(ymin, max(y), ybinsize)
    else:
        raise Exception("Must specify both of either xbinsize/ybinsize or nxbins/nybins")
    counts = n.zeros((nybins, nxbins), dtype='int')
    # Assign each data point (x/y) to a bin.
    for i in range(len(x)):
        xidx = min(int((x[i] - xmin)/xbinsize), nxbins-1)
        yidx = min(int((y[i] - ymin)/ybinsize), nybins-1)
        counts[yidx][xidx] += 1
    # Create 2D x/y arrays, to match 2D counts array.
    xi, yi = n.meshgrid(xbins, ybins)
    return xi, yi, counts

def count_meanminmax(x, y, z, xbinsize=None, ybinsize=None, nxbins=None, nybins=None):
    # Set up grid for contour/density plot.
    xmin = min(x)
    ymin = min(y)
    if (xbinsize!=None) & (ybinsize!=None):
        xbins = n.arange(xmin, max(x), xbinsize)
        ybins = n.arange(ymin, max(y), ybinsize)
        nxbins = xbins.shape[0]
        nybins = ybins.shape[0]
    elif (nxbins!=None) & (nybins!=None):
        xbinsize = (max(x) - xmin)/float(nxbins)
        ybinsize = (max(y) - ymin)/float(nybins)
        xbins = n.arange(xmin, max(x), xbinsize)
        ybins = n.arange(ymin, max(y), ybinsize)
    else:
        raise Exception("Must specify both of either xbinsize/ybinsize or nxbins/nybins")
    counts = n.zeros((nybins, nxbins), dtype='int')
    runsum = n.zeros((nybins, nxbins), dtype='float')
    runmin = n.zeros((nybins, nxbins), dtype='float')
    runmax = n.zeros((nybins, nxbins), dtype='float')
    # Assign each data point (x/y) to a bin.
    for i in range(len(x)):
        xidx = min(int((x[i] - xmin)/xbinsize), nxbins-1)
        yidx = min(int((y[i] - ymin)/ybinsize), nybins-1)
        counts[yidx][xidx] += 1
        runsum[yidx][xidx] += z[i]
        runmin[yidx][xidx] += min(z[i], runmin[yidx][xidx])
        runmax[yidx][xidx] += max(z[i], runmin[yidx][xidx])
    # Create 2D x/y arrays, to match 2D counts array.
    xi, yi = n.meshgrid(xbins, ybins)
    mean = runsum / float(counts)
    return xi, yi, mean, runmin, runmax
