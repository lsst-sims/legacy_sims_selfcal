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
    print 'ackack', stars.size
    names=['x','y','radius']
    types=[float,float,float]

    xyr=np.zeros(stars.size, dtype=zip(names,types))
    # XXX -- should I just make sure everything is radians to start with?
    xyr['x'],xyr['y'] = gnomonic_project_toxy(stars['ra'],stars['dec'], visit['ra'], visit['dec'])
    xyr['r'] = (x**2+y**2)**0.5

    stars =  rfn.merge_arrays([stars, xyr],  flatten=True, usemask=False)

    return stars

def assignPatches(stars, visit, nPatches=16, radiusFoV=1.8):
    """
    Assign PatchIDs to everything.  Assume that stars have already been projected to x,y
    """
    maxval =  gnomonic_project_toxy(0., np.radians(radiusFoV), 0., 0.)

    # This should move all coords to  0 < x < 2.xmaxval
    x = stars['x'] + maxval
    y = stars['y'] + maxval

    nsides = nPatches**0.5
    px = np.floor(x/nsides)
    py = np.floor(y/nsides)

    patchID = np.zeros(stars.size, dtype=zip(['patchID'],[int]))
    patchID += px + py*nsides + visit['visitID']*nPatches
    stars =  rfn.merge_arrays([stars, patchID],  flatten=True, usemask=False)
    return stars
