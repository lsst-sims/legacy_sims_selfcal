from __future__ import print_function
import numpy as np
import scipy.interpolate as interp

"""
A Raft consists of a set of Detectors, each with a thermal model
"""
class FocalplaneThermalModel:
    def __init__(self, focalPlaneDatafile):
        self.rafts = {}
        """For each line in focalPlaneDatafile, add the corresponding raft"""
        fp = open(focalPlaneDatafile, 'r')
        raftList = fp.readlines()
        fp.close()
        for raftSpec in raftList:
            (raftName, raftDataFile, raftOffset) = raftSpec.split()
            raftOffset = float(raftOffset)  
            self.addRaft(raftName, raftDataFile, raftOffset)
                
    def addRaft(self, raftName, raftDataFile, raftOffset):
        print('adding raft %s %s %f\n' % (raftName, raftDataFile, raftOffset))
        self.rafts[raftName] = RaftThermalModel(raftName, raftDataFile, raftOffset)

    def getTemp(self, raftName, detectorName, x, y):
        if (raftName in self.rafts):
            raft = self.rafts[raftName]
            return raft.getTemp(detectorName, x, y)
        else:
            return None
        
            
class RaftThermalModel:
    def __init__(self, raftName, raftDataFile, raftOffset, detectorShape=(4000, 4072)):
        self.raftName = raftName
        self.raftOffset = raftOffset
        self.detectors = {}
        """For each line in raftDataFile, add the corresponding detector"""
        fp = open(raftDataFile, 'r')
        detectorList = fp.readlines()
        fp.close()
        for detectorSpec in detectorList:
            (detectorName, detectorDataFile, offset) = detectorSpec.split()
            detectorOffset = raftOffset + float(offset)
            self.addDetector(detectorName, detectorDataFile, detectorOffset, detectorShape)
    
    def addDetector(self, detectorName, dataFile, offset, detectorShape):
        self.detectors[detectorName] = DetectorThermalModel(detectorShape, dataFile, offset)

    def getTemp(self, detectorName, x, y):
        """
        Given a detectorName and (x,y) coord, return the temp
        """
 
        if (detectorName in self.detectors):
            return self.detectors[detectorName].getTemp(x, y)
        else:
            return None
            
    def getTemps(self, coordList):
        """
        Given a list of detectorNames and (x,y) coords, return the temps
        """
        temps = ()
        for coord in coordList:
            (detectorName, x, y) = coord
            if (detectorName in self.detectors):
                temps.append(self.detectors[detectorName].getTemp(x, y))
            else:
                temps.append(None)
        
        return temps
        
""" 
A DetectorThermalModel is specified by a file with temperature points on a regularly spaced
   grid
"""
class DetectorThermalModel:
    def __init__(self, shape, dataFile, offset):
        self.shape = shape
        self.offset = offset
        
        self.tempData = np.loadtxt(dataFile) + self.offset
        (nx, ny) = self.tempData.shape
        self.x = np.arange(0, shape[0]+1, float(shape[0])/(nx-1))
        self.y = np.arange(0, shape[1]+1, float(shape[1])/(ny-1))
        self.tempFunction = interp.RectBivariateSpline(self.x, self.y, self.tempData)
        
    
    def getTemp(self, x, y):
        return(self.tempFunction(x, y))