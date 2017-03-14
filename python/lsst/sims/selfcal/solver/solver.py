from builtins import zip
from builtins import range
from builtins import object
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsqr
from lsst.sims.selfcal.utils import fastRead

# Modified from /astro/store/scratch/tmp/yoachim/Research/LSST/Parallel_Solver

class lsqrSolver(object):
    """
    Class to read in the output from genCatalog.py and run the self-calibration solver and write the output.

    Might want to expand this so it can be used to fit an arbitrary number of terms?
    """

    def __init__(self, infile=None, patchOut=None, starOut=None, atol=1e-8, btol=1e-8):
        """
        infile:  ASCII file from genCatalog
        patchOut: filename for saving the patch zeropoints
        starOut: filename for saving the star solutions
        atol: tolerance for the solver
        btol: tolerance for the solver
        """

        self.infile = infile
        self.patchOut = patchOut
        self.starOut = starOut
        self.atol=atol
        self.btol=btol

    def run(self):
        self.readData()
        self.cleanData()
        self.solveMatrix()
        self.writeSoln()

    def readData(self):
        names = ['patchID', 'starID', 'mag', 'magErr', 'radius', 'hpID']
        types = [int,int,float,float,float,int]
        self.observations = fastRead(self.infile, dtype=list(zip(names,types)), delimiter=',')

    def cleanData(self):
        """
        Remove observations that can't contribute to a solution.
        Index remaining stars and patches so they are continuous.
        """
        nStart = 1.
        nEnd = 0.
        while nStart != nEnd:
            nStart = self.observations.size
            self.observations.sort(order='starID')
            # Remove observations if the star was only observed once
            good = np.where( (self.observations['starID']-np.roll(self.observations['starID'],1))*
                             (self.observations['starID']-np.roll(self.observations['starID'],-1)) == 0)
            self.observations = self.observations[good]

            # Remove patches with only one star
            self.observations.sort(order='patchID')

            good = np.where( (self.observations['patchID']-np.roll(self.observations['patchID'],1))*
                             (self.observations['patchID']-np.roll(self.observations['patchID'],-1)) == 0)
            self.observations = self.observations[good]
            nEnd = self.observations.size

        self.observations.sort(order='patchID')

        self.Patches = np.unique(self.observations['patchID'])
        nPatches=np.size(self.Patches)
        self.nPatches = nPatches
        self.nPatches = nPatches
        Patches_index=np.arange(nPatches)
        left = np.searchsorted(self.observations['patchID'], self.Patches)
        right = np.searchsorted(self.observations['patchID'], self.Patches, side='right')
        for i in range(np.size(left)): self.observations['patchID'][left[i]:right[i]]=Patches_index[i]

        # Convert starID to continuous running index to keep matrix as small as possible
        self.observations.sort(order='starID')
        self.Stars=np.unique(self.observations['starID'])
        nStars=np.size(self.Stars)
        self.nStars = nStars
        Stars_index = np.arange(1,nStars+1)
        left = np.searchsorted(self.observations['starID'], self.Stars)
        right = np.searchsorted(self.observations['starID'], self.Stars, side='right')
        for i in range(np.size(left)): self.observations['starID'][left[i]:right[i]]=Stars_index[i]



    def solveMatrix(self):
        nObs=np.size(self.observations)
        #construct sparse matrix
        #A = lil_matrix((nPatches+nStars,np.size(observations['PatchID'])))
        row=np.arange(nObs)
        row=np.append(row,row)
        col = np.append(self.observations['patchID'], self.observations['starID']+np.max(self.observations['patchID']))
        #data = np.append(np.ones(nObs),1./observations['magErr'])
        data = 1./self.observations['magErr']
        data = np.append(data,data)
        b = self.observations['mag']/self.observations['magErr'] #maybe do this in place earlier?  then I can just delete parts of observations earlier to save total memory

        #blast away data now that we have the matrix constructed
        del self.observations

        A = coo_matrix( (data,(row,col)), shape = (nObs,self.nPatches+self.nStars))
        A = A.tocsr()
        # solve Ax = b
        self.solution = lsqr(A,b, show=True, atol=self.atol, btol=self.btol)

    def writeSoln(self):

        result = np.empty(self.Patches.size, dtype=list(zip(['patchID','zp'],[int,float])) )
        result['patchID'] = self.Patches
        result['zp'] = self.solution[0][0:self.nPatches]
        np.savez(self.patchOut, result=result)

        result = np.empty(self.Stars.size, dtype=list(zip(['starID','fitMag'],[int,float])) )
        result['starID'] = self.Stars
        result['fitMag'] = self.solution[0][self.nPatches:]
        np.savez(self.starOut, result=result )
