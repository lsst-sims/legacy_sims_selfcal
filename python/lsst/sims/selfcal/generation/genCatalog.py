import numpy as np
import lsst.sims.maf.db as db

def wrapRA(ra):
    """Wraps RA into 0-360 degrees."""
    ra = ra % 360.0
    return ra

def capDec(dec):
    """Terminates declination at +/- 90 degrees."""
    dec = np.where(dec>90, 90, dec)
    dec = np.where(dec<-90, -90, dec)
    return dec
  

def genCatalog(visits, starsDbAddress, offsets=None, lsstFilter='r', raBlockSize=20., decBlockSize=10,
               nPatches=16, radiusFoV=1.8, verbose=True, seed=42):
    """
    Generate a catalog of observed stellar magnitudes.

    visits:  A numpy array with the properties of the visits.  Expected to have Opsim-like values
    starsDbAddress:  a sqlAlchemy address pointing to a database that contains properties of stars used as input.
    offsets:  A list of instatiated classes that will apply offsets to the stars
    lsstFilter:  Which filter to use for the observed stars.
    """

    if offsets is None:
        # Maybe change this to just run with a default SNR offset
        warnings.warn("Warning, no offsets set, returning without running")
        return

    # Set up connection to stars db:
    msrgbDB = db.Database(starsDbAddress, dbTables={'stars':['stars', 'id']})
    starCols = ['id', 'rmag', 'gmag']
    if lsstFilter+'mag' not in starCols:
        starCols.append(lsstFilter+'mag' )

    # Loop over the sky
    raBlocks = np.arange(0.,2.*np.pi, np.radians(raBlockSize))
    decBlocks = np.arange(np.pi, -np.pi, np.radians(decBlockSize) )
    
    for raBlock in raBlocks:
        for decBlock in decBlocks:
            np.random.seed(seed)
            seed += 1 #This should make it possible to run in parallel and maintain repeatability.
            visitsIn = visits(np.where( visits['ra'] >=
                      raBlock & visits['ra'] < raBlock+raBlockSize &
                      visits['dec'] >=  decBlock &
                      visits['dec'] < decBlock+decBlockSize))
            if np.size(visitsIn) > 0:
                # Fetch the stars in this block+the radiusFoV
                decPad = radiusFoV
                raPad = radiusFoV 
                # Need to deal with wrap around effects.
                decMin = capDec(decBlock-decPad)
                decMax = capDec(decBlock+decBlockSize+decPad)
                sqlwhere = 'decl > decMin and decl < decMax '
                raMin = raBlock-raPad/np.cos(np.radians(decMin))
                raMax = raBlock+raBlockSize+raPad/np.cos(np.radians(decMin))

                if wrapRA(raMin) != raMin & wrapRA(raMax) != raMax:
                    # near a pole, just grab all the stars
                    sqlwhere += '' 
                else:
                    raMin = wrapRA(raMin)
                    raMax = wrapRA(raMax)
                    # no wrap
                    if raMin < raMax:
                        slqwhere += 'and ra < raMax and ra > raMin '
                    # One side wrapped
                    else:
                        sqlwhere += 'and (ra > raMin or ra < raMax)'

                stars = msrgbDB.tables['stars'].query_columns_Array(
                    colnames=starCols, constraint=sqlwhere)
                # Add any interesting columns, maybe primary healpix id and hpid 1-4
            for visit in visitsIn:
                # Calc x,y, radius for each star, crop off stars outside the FoV
                # XXX - plan to replace with code to see where each star falls and get chipID, etc
                starsIn = starsProject(starsIn, visit)
                starsIn = starsIn[np.where(starsIn['radius'] <= np.radians(radiusFoV))]

                # Assign patchIDs and healpix IDs
                starsIn = assignPatches(starsIn, visit, nPatches=nPatches)
                #maybe make dmags a dict on stars, so that it's faster to append them all?
                
                # Apply the offsets that have been configured
                for offset in offsets:
                    starsIn = offset.run(starsIn, visit)
                # XXX--print the stars in to a file
                # patchID, starID, observed Mag, mag uncertainty, radius, healpixIDs

                # Note the new starID's and print those to a truth file
                # starID true mag

                # Calc and print a patch file.  Note this is slightly ambiguous since the clouds can have structure
                # patchID magOffset

                # Look at the distribution of dmags
                # patchID, starID, dmag1, dmag2, dmag3...
                
    
