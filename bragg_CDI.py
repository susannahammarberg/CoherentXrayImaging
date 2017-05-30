# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 11:19:26 2017

@author: Susanna Hammarberg
"""
#from IPython import get_ipython
#get_ipython().magic('reset -sf')
#

import numpy as np
import h5py
import matplotlib.pyplot as plt

# Path for raw data from detector lambda
data_path = 'F:\\run48_Wallentin\data\detectors\lambda\\'
sample =  'JWMK16_NW1_cdi1'

nbr_scans = 958#  958

# create matrix to hold diffraction patterns
diffSet=np.zeros((nbr_scans, 516, 1556), dtype=np.int32)

# read data from hdf5-files
for scan_nbr in range(0,nbr_scans): 
    scanCDI = h5py.File(data_path + sample + '\\' + sample + '_' + str('{0:05}'.format(scan_nbr)) + '.nxs','r') # read-only
    data = scanCDI.get('/entry/instrument/detector/data')
    
    np_data = np.array(data)
    
    # rotate and remove 3D thing
    #np_data33 = (np_data33[0])
 #       scan33 = h5py.File('scan33/pilatus_scan_33_' + str('{0:04}'.format(scan_nbr)) + '.hdf5','r') # read-only
 #   data_scan33 = scan33.get('/entry_0000/measurement/Pilatus/data' )
 #   np_data33 = np.array(data_scan33)   #Varför har den tre dimensioner?
    diffSet[scan_nbr] = np_data
 
#description = scanCDI.get(' /entry/instrument/detector/geometry/description')
#np_description = np.array(description)
#del 


#diffSet = misc.imread('P.png',flatten=True)


# Experment parameters
xDet = 0.5849              # Object-detector distance, m; 0.4357 ; 0.5849 ; 0.8919 ; 1.4804
pixel = 0.055E-3           # Pixel ctr to ctr distance (w)
energy = 13.8              # keV
wavelength = 1.23984E-9/energy
kLength = 2*np.pi/wavelength      # 1/m   needed?


plt.figure()
plt.imshow(np.log10(sum(diffSet)))