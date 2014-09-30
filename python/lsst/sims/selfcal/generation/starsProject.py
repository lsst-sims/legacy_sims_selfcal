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
    types=[float]*3
    xyr=np.zeros( dtype=zip(names,types))
    # XXX -- should I just make sure everything is radians to start with?
    xyr['x'],xyr['y'] = gnomonic_project_toxy(stars['ra'],stars['dec'], visit['ra'], visit['dec'])
    xyr['r'] = (x**2+y**2)**0.5

    stars =  rfn.merge_arrays([stars, xyr],  flatten=True, usemask=False)

    return stars
