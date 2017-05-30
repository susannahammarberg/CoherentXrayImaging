# -*- coding: utf-8 -*-
"""
Created on Tue May 30 15:12:51 2017

@author: Sanna
"""
#from IPython import get_ipython
#get_ipython().magic('reset -sf')   #removes all variables saves
import sys   #to collect system path ( to collect function from another directory)
sys.path.insert(0, 'C:/Users/Sanna/Desktop/CXI/Shrinkwrap')

import numpy as np
from numpy import fft
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from create3Dgaussian import create3Dgaussian
from Shrinkwrap3D import shrinkwrap3D 
from shrinkwrap import shrinkwrap 
from math import ceil
from scipy import misc

plt.close("all") # close all plotting windows

# Create 3D crystal: skewed and unskewed.
# create FFTs of these. Filter out one peak.
################################################
Nx = 201
Ny = 201
Nz = 201

crystal3D = np.zeros((Nz,Ny,Nx))
crystal3D_phase = np.zeros((Nz,Ny,Nx))
crystal3D_filter = np.zeros((Nz,Ny,Nx), dtype= np.int32)

dx_filter = 1

#==============================================================================
# for row in range(20,180,dx):  
#     for col in range(80,120,dx):
#         for time in range(90,110,dx): 
#==============================================================================
#def create_3Dcrystal():
index=0
dz = 3
dy = 3
dx = 3

dzrange = [80, 120]
dyrange = [80, 120]
dxrange = [90, 110]

lattice_points_z = np.ceil((dzrange[1] - dzrange[0]) / dz)
lattice_points_y = np.ceil((dyrange[1] - dyrange[0]) / dy)
lattice_points_x = np.ceil((dxrange[1] - dxrange[0]) / dx)             

lattice_points = int(lattice_points_z * lattice_points_y *lattice_points_x)

# coordinates of crystal points (for 3D scattering plot)
#(these are for 3Dscattering plot) makes an index for every 'lattice point'
ycoor = np.zeros((lattice_points))
xcoor = np.zeros((lattice_points))
zcoor = np.zeros((lattice_points))

for depth in range(dzrange[0], dzrange[1], dz):      
    for row in range(dyrange[0], dyrange[1], dy):  
        for col in range(dxrange[0], dxrange[1], dx):
            
            zcoor[index] = depth
            ycoor[index] = row    
            xcoor[index] = col
                        
            crystal3D[depth, row,col,] = 1
            # add phase
            if np.random.randint(2) == 1:
                crystal3D_phase[depth, row, col] = -np.pi/1
            else:
                crystal3D_phase[depth, row, col] = np.pi/1
            index = index+1 


# make an object of both phase and amplitude            
crystal3D = crystal3D * np.exp(1j*crystal3D_phase)
    #return (crystal3D,zcoor,ycoor,xcoor)

def skew_3Dcrystal(crystal, theta):    
    #dx_filter = 1            
    crystal_skewed = np.zeros(( int(np.floor((crystal.shape[0]/np.cos(theta)))) + 50 , crystal.shape[1], crystal.shape[2] ), dtype=np.complex64)


    theta = 10*np.pi/180

    for i in range(0,crystal.shape[0]-6):
        for j in range(0,crystal.shape[1]-6):
            for k in range(0, crystal.shape[2]-6):
                 
                zs = ceil(j / (1 + np.tan(theta)**2) + i*np.tan(theta)/ ( 1 + np.tan(theta)**2)  )
                ys = i
                xs = k
                crystal_skewed[zs,ys,xs] = crystal[i,j,k]
                #np.disp(crystal[i,j,k])
                #if crystal[i,j,k]

    
    crystal_filter2D_skewed = 1
    return (crystal_skewed,  crystal_filter2D_skewed)

theta = 10*np.pi / 180 #rad
crystal3D_skewed, crystal_filter2D_skewed = skew_3Dcrystal(crystal3D, theta )


#diffPattern = abs( fft.fftshift(fft.fftn(crystal_skewed)))**2


#==============================================================================
# plotting
#==============================================================================
def plot_3Dscatter():
#    #scatter 3 ritar ut en ring för varje punkt specificerad av vektorerna (x,y,z)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xcoor,ycoor,zcoor, alpha=0.2, c='red', marker ='.')
#    ax.set_xlim3d(50,150)
#    ax.set_ylim3d(50,150)
#    ax.set_zlim3d(50,150)
plot_3Dscatter()     
                     
 

def plot_crystal3D():
    plt.figure()
    plt.subplot(221)
    plt.imshow(abs(crystal3D[:,:,int(Nx/2)+2]), cmap='gray')
    #plt.title('xyplane cut')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.subplot(222)
    plt.imshow(abs(crystal3D[int(Nz/2)+1,:,:]), cmap='gray')
    #plt.title('xz cut')
    plt.xlabel(' z')
    plt.ylabel(' x')
    plt.colorbar()
    
    plt.subplot(223)        #x
    plt.imshow(abs(crystal3D[:,int(Ny/2)+1,:]), cmap='gray')
    #plt.title('yz cut')
    plt.xlabel(' z')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.suptitle('Plane cuts of crystal amplitude')
        
    plt.figure()   #Phase
    plt.subplot(221)
    plt.imshow(np.angle(crystal3D[:,:,int(Nx/2)+2]), cmap='gray')
    #plt.title('xyplane cut')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.subplot(222)
    plt.imshow(np.angle(crystal3D[int(Nz/2)+1,:,:]), cmap='gray')
    #plt.title('xz cut')
    plt.xlabel(' z')
    plt.ylabel(' x')
    plt.colorbar()
    
    plt.subplot(223)        #x
    plt.imshow(np.angle(crystal3D[:,int(Ny/2)+1,:]), cmap='gray')
    #plt.title('yz cut')
    plt.xlabel(' z')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.suptitle('Plane cuts of crystal phase')
    
plot_crystal3D()    
 
 
def plot_crystal3D():
    plt.figure()
    plt.subplot(221)
    plt.imshow(abs(crystal3D_skewed[:,:,int(Nx/2)+2]), cmap='gray')
    #plt.title('xyplane cut')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.subplot(222)
    plt.imshow(abs(crystal3D_skewed[int(Nz/2)+1,:,:]), cmap='gray')
    #plt.title('xz cut')
    plt.xlabel(' z')
    plt.ylabel(' x')
    plt.colorbar()
    
    plt.subplot(223)        #x
    plt.imshow(abs(crystal3D_skewed[:,int(Ny/2)+1,:]), cmap='gray')
    #plt.title('yz cut')
    plt.xlabel(' z')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.suptitle('Plane cuts of skewed crystal amplitude')
        
    plt.figure()   #Phase
    plt.subplot(221)
    plt.imshow(np.angle(crystal3D_skewed[:,:,int(Ny/2)+2]), cmap='gray')
    #plt.title('xyplane cut')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.subplot(222)
    plt.imshow(np.angle(crystal3D_skewed[int(Ny/2),:,:]), cmap='gray')
    #plt.title('xz cut')
    plt.xlabel(' z')
    plt.ylabel(' x')
    plt.colorbar()
    
    plt.subplot(223)        #x
    plt.imshow(np.angle(crystal3D_skewed[:,int(Ny/2)+1,:]), cmap='gray')
    #plt.title('yz cut')
    plt.xlabel(' z')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.suptitle('Plane cuts of skewed crystal phase')
    
plot_crystal3D()