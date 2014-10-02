import numpy as np
import numpy.lib.recfunctions as rfn



def gnomonic_project_toxy(RA1, Dec1, RAcen, Deccen):
    """Calculate x/y projection of RA1/Dec1 in system with center at RAcen, Deccen.
    Input radians."""
    # also used in Global Telescope Network website
    cosc = np.sin(Deccen) * np.sin(Dec1) + np.cos(Deccen) * np.cos(Dec1) * np.cos(RA1-RAcen)
    x = np.cos(Dec1) * np.sin(RA1-RAcen) / cosc
    y = (np.cos(Deccen)*np.sin(Dec1) - np.sin(Deccen)*np.cos(Dec1)*np.cos(RA1-RAcen)) / cosc
    return x, y

def starsProject(stars, visit):
    """
    Project the stars to x,y plane for a given visit.
    """
    names=['x','y','radius']
    types=[float,float,float]
    stars['x'], stars['y'] = gnomonic_project_toxy(np.radians(stars['ra']),np.radians(stars['decl']),
                                               visit['ra'], visit['dec'])
    stars['radius'] = (stars['x']**2+stars['y']**2)**0.5
    return stars

def assignPatches(stars, visit, nPatches=16, radiusFoV=1.8):
    """
    Assign PatchIDs to everything.  Assume that stars have already been projected to x,y
    """
    maxx, maxy =  gnomonic_project_toxy(0., np.radians(radiusFoV), 0., 0.)
    nsides = nPatches**0.5
    
    # This should move all coords to  0 < x < nsides-1
    px = np.floor((stars['x'] + maxy)/(2.*maxy)*nsides)
    py = np.floor((stars['y'] + maxy)/(2.*maxy)*nsides)
    
    stars['subPatch'] = px + py*nsides 
    stars['patchID'] = stars['subPatch'] + visit['visitID']*nPatches
    return stars
