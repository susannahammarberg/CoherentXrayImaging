# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 15:40:49 2017

@author: Susanna Hammarberg

See s.50 Giewekemeyer thesis för att gå mellan detektorplan och objektplan
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')

# supports sparse matrices. dont know how they work with np matrises
#import scipy.sparse as sparse
from numpy import fft
from ePIE import ePIE  
from create2Dgaussian import create2Dgaussian
import matplotlib.pyplot as plt
import numpy as np
import h5py
from math import floor
from math import ceil
from scipy import misc # to imoport image

from sys import getsizeof   #se hur mkt minne variabler tar upp 
#import matplotlib.animation as animation



plt.close("all")
#directory = 'F:/Nanomax/Wallentin/JWX31C_1/pilatus_scan_49_'
directory = 'F:/Nanomax/Wallentin/JWX31C_1/pilatus_scan_51_'
#directory = 'F:/Nanomax/Wallentin/JWX31C_1/pilatus_scan_49_'
#directory = 'F:/Nanomax/Vogt_ptycho/scan33/pilatus_scan_33_'
#directory = 'F:/Nanomax/Vogt_ptycho/scan83/pilatus_scan_83_'

#metadata_directory = 'F:/Nanomax/Vogt_ptycho/scan83/DiWCr4_2.h5' #2 metadatafiler, den första är 1kB
#metadata_directory = 'F:/Nanomax/Vogt_ptycho/scan33/DiWCr4_1.h5' 
metadata_directory ='F:/Nanomax/Wallentin/JWX31C_1/JWX31C_1.h5' 

# No mask created for Wallentin. (detector is same so mask is the same?) ((No, there are 2 different mask for 2 Voigt scans))
#mask_directory = 'scan33_mask.hdf5'
mask_directory = 'scan83_mask.hdf5'

motorpositions_directory = '/entry51'   
#motorpositions_directory = '/entry49'   
#motorpositions_directory = '/entry83'  
#motorpositions_directory = '/entry33'  


nbr_scans = 440 # 440 Scan49 J.W    #960 #(+1)
nbr_scansy = 21# 31
nbr_scansx = 21#31

# parameters for conversion between detector and object plane
energy = 10.7 # Wallentin 8.17# vogt # keV  
wavelength = (1.23984E-9)/energy
pixel_det = 172E-6   # Pixel ctr to ctr distance (w) [m] #Pilatus
z = 4.3
Ny = 164
Nx = 184#    #192 vogt33   184 Wallentin51
yfactor = (1/Ny)*z*wavelength
xfactor = (1/Nx)*z*wavelength

# create matrix to hold raw diffraction patterns
diffSet=np.zeros((nbr_scans, 195, 487))   

# read data from hdf5-files
for scan_nbr in range(0,nbr_scans): 
    scan = h5py.File( directory  + str('{0:04}'.format(scan_nbr)) + '.hdf5','r') # read-only
    data_scan = scan.get('/entry_0000/measurement/Pilatus/data' )
    np_data = np.array(data_scan)   #Varför har den tre dimensioner?
    # rotate and remove 3D thing
    np_data = (np_data[0])
    diffSet[scan_nbr] = np_data
 

# gather mask

#mask_file = h5py.File(mask_directory)
#meta_mask = mask_file.get('/mask')
#probe_mask = np.array(meta_mask)

def create_mask():
##    probe_mask = np.ones((diffSet.shape[1],diffSet.shape[2]))
    # Find all cold? pixels (they show minus values)
    sumTotalDiffSet= sum(diffSet)
    probe_mask = sumTotalDiffSet > 0
    
    #probe_mask[6,59] = 0
    #probe_mask[6,60] = 0
    #probe_mask[2,243] = 0
    #probe_mask[2,244] = 0
    #probe_mask[11,304] = 0
    #probe_mask[11,305] = 0
    #probe_mask[56,101] = 0
    
    # remove too high intensity pixels
#    probe_mask[111,241] = 0
#    probe_mask[112,240] = 0

    j=239

    probe_mask[111,j:245] = 0
    probe_mask[112,j:245] = 0
    probe_mask[113,j:245] = 0
    probe_mask[114,j:245] = 0
    probe_mask[115,j:245] = 0
    
#    
#    probe_mask[112,241] = 0
#    probe_mask[112,242] = 0
#    probe_mask[112,243] = 0
#    probe_mask[113,241] = 0 
#    probe_mask[113,242] = 0 #pixel with highest intensity
#    probe_mask[113,243] = 0
#    probe_mask[114,241] = 0
#    probe_mask[114,242] = 0
#    probe_mask[114,243] = 0

    return probe_mask

probe_mask = create_mask()

# apply mask
diffSet = probe_mask * diffSet

#del probe_mask

# Trim and center the diffraction patterns around max intensity
# look att the sum of all patterns:
#def trim_center(diffSet):
summed_diffSet = sum(diffSet)
vect = sum(summed_diffSet/summed_diffSet)

### plot sum of diffPatterns
plt.figure()
plt.imshow((sum(diffSet)) , cmap='gray')
plt.title('All diffraction patterns summed')
plt.colorbar()

# inte riktigt rätt för Wallentin, max är vid y=82 (rätt) x=92 (hä, 96)
#diffSet = diffSet[:, 31:195, 150:342] # Vogt33
diffSet = diffSet[:, 31:195, 152:336] # Wallentin51 ändra NxNy också

###plot the trimmed and centered sum of diffPatterns
plt.figure()
plt.imshow(np.log10(sum(diffSet)))
plt.title('All diffraction patterns summed')
plt.colorbar()

# gather motor postions
metadata = h5py.File( metadata_directory)
dataset_motorpositiony = metadata.get(motorpositions_directory + '/measurement/samy')
dataset_motorpositionx = metadata.get(motorpositions_directory + '/measurement/samx')
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
sigmay = 10# 14.1               # initial value of gaussian height
sigmax = 11.8                     # initial value of gaussian width
probe = create2Dgaussian( sigmay, sigmax, diffSet.shape[1], diffSet.shape[2])

def circular_probe():
    circle = misc.imread('circle.png',flatten=True)
    low_values_indices = circle < 200  # Where values are low
    circle[low_values_indices] = 0  # All low values set to 0
    high_values_indices = circle > 0
    circle[high_values_indices] = 1
    return circle
phase = np.pi/2 * np.logical_not(circular_probe())
#phase = np.zeros((diffSet.shape[1],diffSet.shape[2]))
 ##initial guess for probe (square with amplitude 1 everywhere):
#probeInnerSize = 12#with ones
#probe = np.zeros((diffSet.shape[1],diffSet.shape[2]), dtype=np.complex64)
#probe[116:141, 116:141] = 1  #np.ones(shape=(20, 20)). OK
## kolla om detta blev rätt
#probe[floor(diffSet.shape[1]/2 - probeInnerSize/2) :floor(diffSet.shape[1]/2 - probeInnerSize/2)+probeInnerSize,  floor(diffSet.shape[2]/2 - probeInnerSize/2):floor(diffSet.shape[2]/2 - probeInnerSize/2)+probeInnerSize] = 1

#probeInnerSize = 12#with ones
#phase = np.zeros((diffSet.shape[1],diffSet.shape[2]), dtype=np.complex64)
#phase[floor(diffSet.shape[1]/2 - probeInnerSize/2) :floor(diffSet.shape[1]/2 - probeInnerSize/2)+probeInnerSize,  floor(diffSet.shape[2]/2 - probeInnerSize/2):floor(diffSet.shape[2]/2 - probeInnerSize/2)+probeInnerSize] = np.pi/2

## Create square phase 
#phase = np.pi * np.ones((diffSet.shape[1],diffSet.shape[2]), dtype=np.complex64)
# create complex probe?
probe = probe * np.exp(1j*phase)

# size of one pixel in objectplane. (blir annorlunda för att Nx och Ny är olika)
xpixel = xfactor/pixel_det
ypixel = yfactor/pixel_det
# what the width of the diffraction pattern equals to in object plan (pixlar * pixelstorlekx/y)
sizeDiffObjectx =  Nx * xpixel
sizeDiffObjecty =  Ny * ypixel

# hur långt motorn rör sig i x och yled: 
motorWidthx = ( motorpositionx.max() - motorpositionx.min() ) * 1E-6
motorWidthy = ( motorpositiony.max() - motorpositiony.min() ) * 1E-6

# so the size of the object function should be enough to contain:
objectFuncSizeMaxy = motorWidthy + sizeDiffObjecty
objectFuncSizeMaxx = motorWidthx + sizeDiffObjectx

# so with a pixel-size of xpixel * ypixel, the obect function should be this many pixels:
    # should i use ceil!? or floor?
objectFuncNy = ceil(objectFuncSizeMaxy / ypixel)
objectFuncNx = ceil(objectFuncSizeMaxx / xpixel)
# allocate memory for object function
objectFunc = np.zeros((objectFuncNy, objectFuncNx))


# 'normalized' motorpostions converted to meters
positiony = (motorpositiony - motorpositiony.min() ) *1E-6
positionx = (motorpositionx - motorpositionx.min() ) *1E-6


# ska dessa användas ? vet ej men får samma (!?) resultat
#positiony = abs(motorpositiony - motorpositiony.max() ) *1E-6
#positionx = abs(motorpositionx - motorpositionx.max() ) *1E-6

# mirror diffraction patterns
#diffSet = np.fliplr(diffSet)

plt.figure()
plt.imshow(abs(probe), cmap='gray', interpolation='none', extent=[0,3.7943696450428395,0,3.7943696450428395])
plt.xlabel(' [µm]')
plt.ylabel(' [µm]')
plt.title('Probe amplitude')
plt.colorbar()
plt.figure()
plt.imshow(np.angle(probe), cmap='gray', interpolation='none', extent=[0,3.7943696450428395,0,3.7943696450428395])
plt.xlabel(' [µm]')
plt.ylabel(' [µm]')
plt.title('Probe phase')
plt.colorbar()

# run ePIE
objectFunc, probe, ani, sse, psi = ePIE(diffSet, probe, objectFuncNy, objectFuncNx, ypixel, xpixel, positiony, positionx, nbr_scans)

### make ePIE function return psi (exit wave) for every position of probe
# sum over all positions?
psi = sum(psi)
fft1 = fft.fft(psi)
plt.figure()
plt.plot(fft1)

##############################PLOTTING################
plt.show()

#def plott():  
plt.figure()     #, origin="lower"                         # sets the scale on axes. Should be calculated for every new experiment
plt.imshow((np.angle(objectFunc)), interpolation='none', extent=[0,6.837770297837617,0,6.825238081022181])
#plt.gca().invert_yaxis() 
plt.xlabel(' [µm]')
plt.ylabel(' [µm]')
plt.title('Object phase')
#plt.clim(-np.pi,np.pi)
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
#plt.clim(-np.pi,np.pi)
plt.colorbar()

plt.figure()
plt.plot(sse)
plt.xlabel(' iterations ')
plt.ylabel(' SSE ')
plt.title('SSE')

plt.figure()
plt.plot(abs(probe.sum(axis=0)))
plt.xlabel(' [µm]')
plt.title('Probe summed over all columns')
#plt.axis((0,xpixel*diffSet.shape[2]*1E6,0,10))

#my_xticks = xpixel*diffSet.shape[2]*1E6
#plt.set_xticks( my_xticks )

plt.figure()
plt.plot(abs(probe.sum(axis=1)))  #sum 
plt.xlabel(' [µm]')
plt.title('Probe summed over all rows')

#plt.axis((0,ypixel*diffSet.shape[1]*1E6,0,20))
 #   return 0


#def get_size(obj, seen=None):
#    """Recursively finds size of objects"""
#    size = getsizeof(obj)
#    if seen is None:
#        seen = set()
#    obj_id = id(obj)
#    if obj_id in seen:
#        return 0
#    # Important mark as seen *before* entering recursion to gracefully handle
#    # self-referential objects
#    seen.add(obj_id)
#    if isinstance(obj, dict):
#        size += sum([get_size(v, seen) for v in obj.values()])
#        size += sum([get_size(k, seen) for k in obj.keys()])
#    elif hasattr(obj, '__dict__'):
#        size += get_size(obj.__dict__, seen)
#    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
#        size += sum([get_size(i, seen) for i in obj])
#    return size
#