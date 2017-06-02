# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 13:20:16 2017

@author: Sanna
"""
import numpy as np
from numpy import fft
from math import ceil
import matplotlib.pyplot as plt

# Nbr of used pixels on the detector
Nxy_det = 190
energy = 10.72   #keV    
wavelength = 1.23984E-9 / energy
pixel_det = 172E-6   # Pixel ctr to ctr distance (w) [m] #Pilatus
z = 5


def calculate_dq():
    # lattice constant InP
    lattice_constant_a = 5.8687E-10 
    # distance between planes (111), correct?
    d = lattice_constant_a / np.sqrt(3) 
    q_abs = 2*np.pi / d
    
    dtheta = 0.001 #degrees
    dq = q_abs * dtheta

# one pixel on the object plane correspongs to ... m in the object plane at distance z from the detector
pixel_objectPlane = wavelength * z / ( Nxy_det * pixel_det)

# the field of view in the object plane is
FOV_real_size = Nxy_det * pixel_objectPlane

# FOV in the object plane at distance z (with pixel sizes pixel_objectPlane)
FOV_l = np.zeros((Nxy_det,Nxy_det))

# size of real object in [m]
x_obj = 1000E-9
y_obj = 100E-9
# how many pixels (object plane pixels) that the object sizes fits in 
Nx_obj = ceil( x_obj / pixel_objectPlane)
Ny_obj = ceil( y_obj / pixel_objectPlane)

# create object of size N_obj 
#obj = np.ones((Ny_obj,Nx_obj))
#or a 2D crystal object:
obj = np.zeros((Ny_obj,Nx_obj))

dy = 2
dx = 2

dyrange = [0, Ny_obj]
dxrange = [0, Nx_obj]

for row in range(dyrange[0], dyrange[1], dy):  
    for col in range(dxrange[0], dxrange[1], dx):
        obj[row,col] = 1

# insert the object into the FOV
FOV_l[int((Nxy_det-Ny_obj)/2) : int((Nxy_det-Ny_obj)/2) + Ny_obj, int((Nxy_det-Nx_obj)/2) : int((Nxy_det-Nx_obj)/2) + Nx_obj ] = obj
      
# FFT2 of the space in FOV at that distance z      
# this is what the pattern from the obj would look lika at that distance
fft_FOV_l = fft.fftshift(fft.fft2(FOV_l))      

plt.figure()
plt.subplot(221)
pixel_objectPlane_nano = pixel_objectPlane *1E9
plt.title('z=0, dx=%.2f nm '%pixel_objectPlane_nano)
plt.imshow(FOV_l, extent=[0, FOV_real_size*1E6, 0, FOV_real_size*1E6])#, cmap='gray')  FOV_real_size
plt.xlabel(' x [um]')
plt.ylabel(' y [um]')


plt.subplot(222)
plt.title('z=%.2f m dx = 172 um'%z)
plt.imshow(abs(fft_FOV_l), extent=[0, Nxy_det * pixel_det*1E3, 0, Nxy_det * pixel_det*1E3])#, cmap='gray')  FOV_real_size
plt.xlabel(' x [mm]')
plt.ylabel(' y [mm]')

plt.savefig('savefig\diff_sim_crystal_%d'%z, bbox_inches='tight')
