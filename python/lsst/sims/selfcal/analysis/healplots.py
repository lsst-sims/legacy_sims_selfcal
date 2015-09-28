import numpy as np
import healpy as hp

def healbin(ra, dec, values, nside=128, reduceFunc=np.mean, dtype='float', returnPixAng=False):
    """Take arrays of ra, dec, and value and bin into healpixels  """

    lat = np.pi/2. - dec
    hpids = hp.ang2pix(nside, lat, ra)

    order = np.argsort(hpids)
    hpids = hpids[order]
    values = values[order]
    pixids = np.arange(hp.nside2npix(nside))

    left = np.searchsorted(hpids, pixids)
    right = np.searchsorted(hpids, pixids, side='right')

    mapVals = np.zeros(pixids.size, dtype=dtype)

    for idx in pixids:
        mapVals[idx] = reduceFunc(values[left[idx]:right[idx]] )


    if returnPixAng:
        lat, lon = hp.pix2ang(nside,pixids)
        decRet = np.pi/2. - lat
        raRet = lon
        return mapVals, raRet, decRet
    else:
        return mapVals
