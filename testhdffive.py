# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 15:40:49 2017

@author: HonkyT
"""
# supports sparse matrices. dont know how they work with np matrises
import scipy.sparse as sparse
from ePIE import ePIE  
import matplotlib.pyplot as plt
import numpy as np
import h5py
from math import floor

from sys import getsizeof   #se hur mkt minne variabler tar upp 

nbr_scans = 960
nbr_scansx = 31
nbr_scansy = 31
detectorSizey = 195 #487
detectorSizex = 487 # 195

#conversionfactors between detector plane and object plane
# one for xled one for yled. 
# xfactor = lambda*z / Nx
wavelength = (1E-9)*1.24/8170
z = 4.3
Nx = 200
Ny = 195
xfactor = (1/Nx)*z*wavelength
yfactor = (1/Ny)*z*wavelength
# allocate memory for the diffraction patterns
diffSet=np.zeros((nbr_scans, detectorSizey, detectorSizex)) #have to uese double paranthesis   , dtype=complex

# get size of the variabel
#getsizeof(diffSet)

for scan_nbr in range(0,nbr_scans): 
    scan33 = h5py.File('scan33/pilatus_scan_33_' + str('{0:04}'.format(scan_nbr)) + '.hdf5','r') # read-only
    data_scan33 = scan33.get('/entry_0000/measurement/Pilatus/data' )
    np_data33 = np.array(data_scan33)   #Varför har den tre dimensioner?
    # rotate and remove 3D thing
    np_data33 = (np_data33[0])
    
    diffSet[scan_nbr] = np_data33
    
    #plt.figure()
    #plt.imshow(np.log10(diffSet[scan_nbr]))

addOrigDiffSet = np.log10(sum(diffSet))

# Trim the diffraction patterns in x-led
#diffSet = np.zeros((nbr_scans, 200, ))
diffSet = diffSet[:,:,150:350]
# to looka at the hole diffSet
addDiffSet = np.log10(sum(diffSet))
# gather motor postions
metadata = h5py.File('DiWCr4_1.h5')
dataset_motorpositiony = metadata.get('/entry33/measurement/samy')
dataset_motorpositionx = metadata.get('/entry33/measurement/samx')
motorpositiony = np.array(dataset_motorpositiony) 
motorpositionx = np.array(dataset_motorpositionx) 

# calculate how long each step is in x and y  OBS kan också vara minus
stepSizex = np.zeros((nbr_scansx,1))
stepSizey = np.zeros((nbr_scansy,1))
for i in range(0,nbr_scansx):   #gör 2 loops for diffrent nbr of scans in y and x 
    stepSizex[i] = motorpositionx[i+1] - motorpositionx[i]    # obs obs detta är i micrometerOBS
    stepSizey[i] = motorpositiony[i+1] - motorpositiony[i]

# gather mask
mask_file = h5py.File('scan33_mask.hdf5')
meta_mask = mask_file.get('/mask')
probe_mask = np.array(meta_mask)
 

# initial guess for probe (the inner size with ones) and outer size as diffSet:
probeInnerSize = 25 #with ones
probe = np.zeros((diffSet.shape[1],diffSet.shape[2]), dtype=np.complex64)
#probe[116:141, 116:141] = 1  #np.ones(shape=(20, 20)). OK
# kolla om detta blev rätt
probe[floor(diffSet.shape[1]/2 - probeInnerSize/2) :floor(diffSet.shape[1]/2 - probeInnerSize/2)+probeInnerSize,  floor(diffSet.shape[2]/2 - probeInnerSize/2):floor(diffSet.shape[2]/2 - probeInnerSize/2)+probeInnerSize] = 1

objectSizey = 500
objectSizex = 500

# run ePIE
animation = ePIE(diffSet, probe, objectSizey, objectSizex)


#test = h5py.File('test_scan.hdf5','r') # read-only  # skickar nu tillbaka objectfunc
#
#data = test.get('/entry_0000/measurement/Pilatus/data' )
#
#np_data=np.array(data)   #Varför har den tre dimensioner?
#
#plt.figure()
#
#np_data_int=np.dtype(float)
#plt.imshow(np.log10(np_data[0])) # du kan itne göra log10 på int32 data


## måste jag skapa grupppen
#grupp = test.create_gru
##f = h5py.File("mytestfile.hdf5", "w")
##dset = f.create_dataset("mydataset", (100,), dtype='i')
##
##grp = f.create_group("subgroup")
##
##dset2 = grp.create_dataset("another_dataset", (50,), dtype='f')
##
##
##dset3 = f.create_dataset('/entry_0000/measurement/Pilatus/data', dtype='i')
##
##
##list_of_item_names = test.items()
##print( list_of_item_names)
###test.name
###test.filenam
###test.keys
##vad=isinstance(test, h5py.File) #or Group or Dataset
##
#
##dset = l.create_dataset("mydataset", (100,), dtype='i')
#
##filenames pilatus_scan_33_0000.hdf5
#
##/* Open an existing dataset. */
#
##   
#for name in entry_0000:
#    print( name)
#    
#    #"entry_0000" in test    ==true!!
#    
#    
#    
##fileName = "prj_test.nexus.hdf5"
#fileName = "test_scan.hdf5"
#f = h5py.File(fileName,  "r")
#
#for item in f.attrs.keys():
#    print( item + ":", f.attrs[item])
#    
##mr = f['/entry/mr_scan/mr']
#mr = f['/entry_0000']
##i00 = f['/entry/mr_scan/I00']
##print( "%s\t%s\t%s" % ("#", "mr", "I00"))
##for i in range(len(mr)):
##    print ("%d\t%g\t%d" % (i, mr[i], i00[i]))
##f.close()
