#
#   Module mysort.py
#
# Purpose : sort arrays by specified column using argsort
#    for multicolumn sorting see lexsort
#
import numpy
import pylab
import random
#      x= numpy.random.randint(-10,10,len(am))

def csort((arr,icol)):
        sar = numpy.zeros(numpy.shape(arr))
        idx = numpy.argsort(arr[icol])
        for i in idx :
                sar[:,i] = arr[:,idx[i]]
        return sar
        
       

# test sorting
#
#  --create array ta
tl = [(1.,2.,3.,4.),(10.,30.,40.,20.),(400.,300.,200.,100.)]
ta = numpy.array(tl)
