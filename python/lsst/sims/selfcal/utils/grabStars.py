import numpy as np

# Utility to grab a fairly uniform density of stars from the fatboy database

# Note, solid angle of Latitude-longitude rectangle = (sin (dec_north) - sin(dec_south)(RA_east - RA_west)


# Number of total stars wanted on whole sphere:
Nstars = 2e6

# Degrees?
decSize = 1.
raSize = 1.

decBlocks = np.arange(-90.,90.,decSize)
raBlocks = np.arange(0.,360., raSize)

cols = ['ra','decl', 'id', 'fe/h', 'logg', 'umag','gmag','rmag',
        'imag','zmag', 'ymag']



for ra in raBlocks:
    for dec in decBlocks:
        solidAngle = (np.sin(np.radians(dec+decSize))-np.sin(np.radians(dec)))*(np.radians(raSize))

        starsInBlock = Nstars * solidAngle / (4.*np.pi)
        
        sqlw = 'where (ra between %f and %f) and (decl between %f and %f)'%(ra, ra+raSize, dec, dec+decSize)

        
