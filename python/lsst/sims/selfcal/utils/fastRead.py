from __future__ import print_function
import numpy as np
import re

def fastRead(infilename, dtype,
             delimiter=None, verbose=True, comments=['!','#','%']):
    """
    Read values from infilename, and return numpy array.
    The numpy reading tools (genfromtxt, loadtxt) are incredibly slow.

    infilename: File to read
    dtype: numpy dtype for the returned structured array
    delimiter: How values are seperated
    comments: If a like starts with one of these characters, skip it.
    """
    # Open file.
    f = open(infilename, 'r')
    if verbose:
        print("# Reading file %s" %(infilename))
    # Read data from file.
    value = []
    for line in f:
        # Ignore line if it's a comment.
        if line[0] in comments:
            continue
        else:
            linevalues = line.strip().split(delimiter)
        value.append(tuple(linevalues))
    f.close()
    # Convert to numpy array.
    value = np.array( value, dtype=dtype)
    return value
