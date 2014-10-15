import numpy as np
import readDatafile as rDf
import scipy as sp
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsqr

#read in data

observations = rDf.readDatafile('star_obs.dat', \
                   ['PatchID', 'StarID', 'mag', \
                    'magErr', 'FpPatchID'], \
                    dtype=[('PatchID', '<i8'), ('StarID', '<i8'), \
                    ('mag','<f8'), ('magErr','<f8'), \
                    ('FpPatchID','<i8')],skiprows=1)

#remove observations if star only observed once

observations = np.sort(observations, order='StarID')
good = np.where( (observations['StarID']-np.roll(observations['StarID'],1))*
                 (observations['StarID']-np.roll(observations['StarID'],-1)) == 0)
observations = observations[good]
Stars=np.unique(observations['StarID'])
nStars=np.size(Stars)
Stars_index = np.arange(nStars)
left = np.searchsorted(observations['StarID'], Stars)
right = np.searchsorted(observations['StarID'], Stars, side='right')
for i in range(np.size(left)): observations['StarID'][left[i]:right[i]]=Stars_index[i]



#remove patches where there is only one star

observations = np.sort(observations, order='PatchID')
good = np.where( (observations['PatchID']-np.roll(observations['PatchID'],1))*
                 (observations['PatchID']-np.roll(observations['PatchID'],-1)) == 0)
observations = observations[good]
Patches = np.unique(observations['PatchID'])
nPatches=np.size(Patches)
Patches_index=np.arange(nPatches)
#convert PatchID to patch index

left = np.searchsorted(observations['PatchID'], Patches)
right = np.searchsorted(observations['PatchID'], Patches, side='right')
for i in range(np.size(left)): observations['PatchID'][left[i]:right[i]]=Patches_index[i]

nObs=np.size(observations)
#construct sparse matrix
#A = lil_matrix((nPatches+nStars,np.size(observations['PatchID'])))
row=np.arange(nObs)
row=np.append(row,row)
col = np.append(observations['PatchID'], observations['StarID']+np.max(observations['PatchID']))
#data = np.append(np.ones(nObs),1./observations['magErr'])
data = 1./observations['magErr']
data = np.append(data,data)
b = observations['mag']/observations['magErr'] #maybe do this in place earlier?  then I can just delete parts of observations earlier to save total memory

#blast away data now that we have the matrix constructed
del observations

A = coo_matrix( (data,(row,col)), shape = (nObs,nPatches+nStars))
A = A.tocsr()

solution = lsqr(A,b,iter_lim=2668, show=True)

#output the bestfit parameters
rDf.writeDatafile('test_bestfit_Patch.dat', {'Patches':Patches, 'ZP':solution[0][0:nPatches]}, ['Patches','ZP'], printheader=False)

rDf.writeDatafile('test_bestfit_Star.dat', {'Stars':Stars, 'Mag':solution[0][nPatches:]}, ['Stars','Mag'], printheader=False)
