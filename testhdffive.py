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
from math import ceil

from sys import getsizeof   #se hur mkt minne variabler tar upp 
import matplotlib.animation as animation

nbr_scans = 960
nbr_scansy = 31
nbr_scansx = 31

# parameters for conversion between detector and object plane
energy = 8.17 # keV
wavelength = (1E-9)*1.24/energy
z = 4.3
Ny = 192
Nx = 164
yfactor = (1/Ny)*z*wavelength
xfactor = (1/Nx)*z*wavelength

# allocate memory for the diffraction patterns
diffSet=np.zeros((nbr_scans, 195, 487))

# get size of the variabel
#getsizeof(diffSet)

# läs in data från hdf5-filerna
for scan_nbr in range(0,nbr_scans): 
    scan33 = h5py.File('scan33/pilatus_scan_33_' + str('{0:04}'.format(scan_nbr)) + '.hdf5','r') # read-only
    data_scan33 = scan33.get('/entry_0000/measurement/Pilatus/data' )
    np_data33 = np.array(data_scan33)   #Varför har den tre dimensioner?
    # rotate and remove 3D thing
    np_data33 = (np_data33[0])
    
    diffSet[scan_nbr] = np_data33
    
# bara för att titta på summan av alla diffSet, originalstorlek
addOrigDiffSet = np.log10(sum(diffSet))

# Trim and cetner the diffraction patterns around max intensity
#diffSet = np.zeros((nbr_scans, 200, ))
diffSet = diffSet[:, 31:195, 150:342]
# bara för att titta på summan av alla diffSet, trimmade
addDiffSet = np.log10(sum(diffSet))
# test-trimmning på summan av diffmönster
#addDiffSet = addDiffSet[31:195, 0:192 ]

# gather motor postions
metadata = h5py.File('DiWCr4_1.h5')
dataset_motorpositiony = metadata.get('/entry33/measurement/samy')
dataset_motorpositionx = metadata.get('/entry33/measurement/samx')
motorpositiony = np.array(dataset_motorpositiony) 
motorpositionx = np.array(dataset_motorpositionx) 

# calculate how long each step is in x and y OBS kan också vara minus
stepSizex = np.zeros((nbr_scansx,1))
stepSizey = np.zeros((nbr_scansy,1))
for i in range(0,nbr_scansx):   #gör 2 loops for diffrent nbr of scans in y and x . convert from microns to meters
    stepSizex[i] = (motorpositionx[i+1] - motorpositionx[i]) * 1E-6
    stepSizey[i] = (motorpositiony[i+1] - motorpositiony[i]) * 1E-6

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

objectSizey = 353
objectSizex = 358

# #test
# size of one pixel in objectplane. (blir annorlunda för att Nx och Ny är olika)
xpixel = xfactor/(172E-6)
ypixel = yfactor/(172E-6)
# what the width of the diffraction pattern equals to in object plan (pixlar * pixelstorlekx/y)
sizeDiffObjectx =  Nx * xpixel
sizeDiffObjecty =  Ny * ypixel

# hur långt motorn rör sig i x och yled: 
motorWidthx = ( motorpositionx.max() - motorpositionx.min()    ) * 1E-6
motorWidthy = ( motorpositiony.max() - motorpositiony.min()  )* 1E-6

# so the size of the object function should be enough to contain:
objectFuncSizeMaxy = motorWidthy + sizeDiffObjecty
objectFuncSizeMaxx = motorWidthx + sizeDiffObjectx

# so with a pixel-size of xpixel * ypixel, the obect function should be this many pixels:
    # should i use ceil!? or floor?
objectFuncNy = ceil(objectFuncSizeMaxy / ypixel)
objectFuncNx = ceil(objectFuncSizeMaxx / xpixel)

objectFunc = np.zeros((objectFuncNy, objectFuncNx))

#for matrixIndex in range(0, 960+1):           #egentligen onödigt med 2 loopar. lättare att läsa
#    
# ska ju matcha med diffraktionsmöstreran, så kolla det

# för matrixIndex = 0:
#objectFunc[ 0: 195 , 0: 200] = 5
# för matrixIndex = 1:
#stepSizey # hur många pixlar motsvarar detta?
# ( eg. oldindexy + np.round(stepSizey[0]/ypixel))
#objectFunc[ int(np.round(stepSizey[0]/ypixel)) : int(np.round(stepSizey[0]/ypixel)) + 195, int(np.round(stepSizex[0]/xpixel)) : int(np.round(stepSizex[0]/xpixel)) + 200 ] = 1#diffSet[0]


fig = plt.figure()

# Initialize vector for animation data
ims = []
testy = (motorpositiony - motorpositiony.min() ) *1E-6
testx = (motorpositionx - motorpositionx.min() ) *1E-6

for u in range(0,961):

    yposition = int(np.round(testy[u]/ypixel))    
    xposition = int(np.round(testx[u]/xpixel))
    objectFunc[yposition : yposition +  Ny, xposition : xposition + Nx ] = objectFunc[0 : Ny, xposition : xposition + Nx] + 1
     # anim
    im = plt.imshow(np.log10(objectFunc), animated=True)
    ims.append([im])
    
    
ani = animation.ArtistAnimation(fig, ims, interval=150, blit=True,repeat_delay=2000)
#aniT = animation.ArtistAnimation(figT, imsT, interval=4000, blit=True,repeat_delay=2000)
plt.show()
    
#plt.figure()
#plt.imshow(np.log10(objectFunc))

     # Cut out the part of the image that is illuminated at R(=(ypos,xpos)
#     objectIlluminated[0] = objectFunc[0: 3.7 um, 0 : 3.8 um]
#                                             med omvanlidingsfaktor
#     objectIlluminated[1] = objectFunc[0 + stepssize[0] : ans + 3.7, 
#                                        (where stepsize = motorposition i+1 - 1 * omvandlingsfaktor)
#     objectIlluminated = objectFunc[ypos*stepsize:ypos*stepsize+ysize, xpos*stepsize:xpos*stepsize+xsize]

# run ePIE
#animation = ePIE(diffSet, probe, objectSizey, objectSizex)


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
