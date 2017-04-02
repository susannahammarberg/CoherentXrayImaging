# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 15:40:49 2017

@author: HonkyT
"""
# supports sparse matrices. dont know how they work with np matrises
#import scipy.sparse as sparse
from ePIE import ePIE  
from create2Dgaussian import create2Dgaussian
import matplotlib.pyplot as plt
import numpy as np
import h5py
#from math import floor
from math import ceil

from sys import getsizeof   #se hur mkt minne variabler tar upp 
#import matplotlib.animation as animation


plt.close("all")


nbr_scans = 960
nbr_scansy = 31
nbr_scansx = 31

# parameters for conversion between detector and object plane
energy = 8.17 # keV
wavelength = (1E-9)*1.24/energy
z = 4.3
Ny = 164
Nx = 192
yfactor = (1/Ny)*z*wavelength
xfactor = (1/Nx)*z*wavelength

# create matrix to hold diffraction patterns
diffSet=np.zeros((nbr_scans, 195, 487))

# read data from hdf5-files
for scan_nbr in range(0,nbr_scans): 
    scan33 = h5py.File('scan33/pilatus_scan_33_' + str('{0:04}'.format(scan_nbr)) + '.hdf5','r') # read-only
    data_scan33 = scan33.get('/entry_0000/measurement/Pilatus/data' )
    np_data33 = np.array(data_scan33)   #Varför har den tre dimensioner?
    # rotate and remove 3D thing
    np_data33 = (np_data33[0])
    
    diffSet[scan_nbr] = np_data33
 

# gather mask
mask_file = h5py.File('scan33_mask.hdf5')
meta_mask = mask_file.get('/mask')
probe_mask = np.array(meta_mask)

# apply mask
diffSet = (probe_mask * diffSet)
del probe_mask

# Trim and center the diffraction patterns around max intensity
diffSet = diffSet[:, 31:195, 150:342]


# gather motor postions
metadata = h5py.File('DiWCr4_1.h5')
dataset_motorpositiony = metadata.get('/entry33/measurement/samy')
dataset_motorpositionx = metadata.get('/entry33/measurement/samx')
motorpositiony = np.array(dataset_motorpositiony) 
motorpositionx = np.array(dataset_motorpositionx) 

del metadata
 
# calculate how long each step is in x and y OBS kan också vara minus
stepSizex = np.zeros((nbr_scansx,1))
stepSizey = np.zeros((nbr_scansy,1))
for i in range(0,nbr_scansx):   #gör 2 loops for diffrent nbr of scans in y and x . convert from microns to meters
    stepSizex[i] = (motorpositionx[i+1] - motorpositionx[i]) * 1E-6
    stepSizey[i] = (motorpositiony[i+1] - motorpositiony[i]) * 1E-6

# probe construction
sigmay = 15.1;               # initial value of gaussian height
sigmax = 11#11.8;               # initial value of gaussian width
probe = create2Dgaussian( sigmay, sigmax, diffSet.shape[1], diffSet.shape[2])

# initial guess for probe (the inner size with ones) and outer size as diffSet:
#probeInnerSize = 12#with ones
#probe = np.zeros((diffSet.shape[1],diffSet.shape[2]), dtype=np.complex64)
#probe[116:141, 116:141] = 1  #np.ones(shape=(20, 20)). OK
# kolla om detta blev rätt
#probe[floor(diffSet.shape[1]/2 - probeInnerSize/2) :floor(diffSet.shape[1]/2 - probeInnerSize/2)+probeInnerSize,  floor(diffSet.shape[2]/2 - probeInnerSize/2):floor(diffSet.shape[2]/2 - probeInnerSize/2)+probeInnerSize] = 1

# size of one pixel in objectplane. (blir annorlunda för att Nx och Ny är olika)
xpixel = xfactor/(172E-6)
ypixel = yfactor/(172E-6)
# what the width of the diffraction pattern equals to in object plan (pixlar * pixelstorlekx/y)
sizeDiffObjectx =  Nx * xpixel
sizeDiffObjecty =  Ny * ypixel

# hur långt motorn rör sig i x och yled: 
motorWidthx = ( motorpositionx.max() - motorpositionx.min()  ) * 1E-6
motorWidthy = ( motorpositiony.max() - motorpositiony.min()  )* 1E-6

# so the size of the object function should be enough to contain:
objectFuncSizeMaxy = motorWidthy + sizeDiffObjecty
objectFuncSizeMaxx = motorWidthx + sizeDiffObjectx

# so with a pixel-size of xpixel * ypixel, the obect function should be this many pixels:
    # should i use ceil!? or floor?
objectFuncNy = ceil(objectFuncSizeMaxy / ypixel)
objectFuncNx = ceil(objectFuncSizeMaxx / xpixel)
# allocate memory for object function
objectFunc = np.zeros((objectFuncNy, objectFuncNx))


#
positiony = (motorpositiony - motorpositiony.min() ) *1E-6
positionx = (motorpositionx - motorpositionx.min() ) *1E-6

# run ePIE
objectFunc, probe, ani = ePIE(diffSet, probe, objectFuncNy, objectFuncNx, ypixel, xpixel, positiony, positionx)

##############################PLOTTING################
plt.show()

#def plott():  
plt.figure()     #, origin="lower"
plt.imshow((np.angle(objectFunc)), interpolation='none', extent=[0,6.837770297837617,0,6.825238081022181])
   # plt.gca().invert_yaxis()  Detta är väl som Alex eg
plt.xlabel(' [µm]')
plt.ylabel(' [µm]')
plt.title('Object phase')
plt.colorbar()
   
plt.figure()                                                            # horisontalt vertikalt. xpixel * size(objectfunc[xled])
plt.imshow(abs(objectFunc), interpolation='none', extent=[0,6.837770297837617,0,6.825238081022181])
plt.xlabel(' [µm]')
plt.ylabel(' [µm]')
plt.title('Object amplitude')
plt.colorbar()
    
plt.figure()
plt.imshow(abs(probe), interpolation='none', extent=[0,3.7943696450428395,0,3.7943696450428395])
plt.xlabel(' [µm]')
plt.ylabel(' [µm]')
plt.title('Probe amplitude')
plt.colorbar()

plt.figure()                                                            # horisontalt vertikalt
plt.imshow(np.angle(probe), interpolation='none', extent=[0,3.7943696450428395,0,3.7943696450428395])
plt.xlabel(' [µm]')
plt.ylabel(' [µm]')
plt.title('Probe phase')
plt.colorbar()
 #   return 0



def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


