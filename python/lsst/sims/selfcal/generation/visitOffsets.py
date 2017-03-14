from builtins import zip
import numpy as np
import numpy.lib.recfunctions as rfn


def visitOffsets(visits, zpOff=1.):
    """
    Convert an opsim array to have more generic names and
    add any extra columns we may want
    """

    dnames = visits.dtype.names
    dnames = ['ra' if (x=='fieldRA' or x=='ditheredRA') else x for x in dnames]
    dnames = ['dec' if (x=='fieldDec' or x=='ditheredDec') else x for x in dnames]
    dnames = ['visitID' if (x=='obsHistID') else x for x in dnames]
    
    visits.dtype.names = dnames
    
    zp = (np.random.uniform(size=visits.size)*zpOff).astype(list(zip(['zpOff'],[float])))
    visits = rfn.merge_arrays([visits, zp],  flatten=True, usemask=False)
    return visits
