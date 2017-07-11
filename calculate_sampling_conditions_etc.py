# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 13:20:16 2017

@author: Sanna
"""
import numpy as np
from numpy import fft
from math import ceil
from scipy.interpolate import RegularGridInterpolator    #to interpolate 3d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

plt.close("all") # close all plotting windows


#-------------------------------------
# This version of the code is used to calculate
# pixel sizes, theta, dthea etc on NanoMAX
#--------------------------------------

# Nbr of used pixels on the detector
N1 = 512 #190   --->z
N2 = 512#25 #190   --->y
N3 = 51
   #30? different angles are used --->x
energy = 9.5   #10.72   #keV    
wavelength = 1.23984E-9 / energy
pixel_det = 55E-6   # Pixel ctr to ctr distance (w) [m] #Pilatus
z = 1

#def calculate_dq3():
# lattice constant Ga(0.51In(0.49)P
#lattice_constant_a = 5.653E-10 
# lattice constant InP
lattice_constant_a = 5.8687E-10 
# distance between planes (111), correct?
d = lattice_constant_a / np.sqrt(3) 
q_abs = 2*np.pi / d
#q_abs_alt = 4*np.pi *np.sin(theta)/wavelength
theta = np.arcsin(q_abs*wavelength/(4*np.pi))
theta = theta * 180/ np.pi
np.disp('Theta [o]:')
np.disp(theta)

dtheta = wavelength/ (4*(500E-9)*np.sin(theta *np.pi/180))
dtheta = dtheta * 180 / np.pi
print 'dtheta [o]' 
np.disp(dtheta)
# reconvert to radians for later calcultations
theta = theta * np.pi / 180
dtheta = dtheta * np.pi / 180
#    return dq

# define pixel sizes in reciprocal space
dq1 = 2*np.pi*pixel_det /(wavelength * z)
dq2 = 2*np.pi*pixel_det /(wavelength * z)
dq3 = q_abs * dtheta
                
q1_linspace = np.linspace( - dq1*N1 /2, dq1*N1 /2, N1 )
q2_linspace = np.linspace( - dq2*N2 /2, dq2*N2 /2, N2 )
q3_linspace = np.linspace( - dq3*N3 /2, dq3*N3 /2, N3 )

q1, q2, q3 = np.meshgrid( q1_linspace, q2_linspace, q3_linspace )

# one pixel on the object plane correspongs to ... m in the object plane at distance z from the detector (x and y)
dr1 = 2*np.pi / ( N1 * dq1 * np.cos(theta))
dr2 = 2*np.pi / ( N2 * dq2 * np.cos(theta))
dr3 = 2*np.pi / ( N3 * dq3 * np.cos(theta))

#pixel_objectPlane = wavelength * z / ( Nxy_det * pixel_det)
#pixel_objectPlane_z = wavelength * z / ( Nz_det * dq_3)

# not used: but good to check with
# the field of view in the object plane is (this is just a number (size of r-space in x and y z direction))
FOV_size_r1 = N1 * dr1
FOV_size_r2 = N2 * dr2
FOV_size_r3 = N3 * dr3

# FOV in the object plane at distance z (with pixel sizes pixel_objectPlane) (this represents the actual space)
FOV_r = np.zeros((N3,N2,N1))

r1_linspace = np.linspace( - dr1*N1 /2, dr1*N1 /2, N1 )
r2_linspace = np.linspace( - dr2*N2 /2, dr2*N2 /2, N2 )
r3_linspace = np.linspace( - dr3*N3 /2, dr3*N3 /2, N3 )

r1, r2, r3 = np.meshgrid(r1_linspace, r2_linspace, r3_linspace)

# create realspace pixel sizes
dx = dr3*np.cos(theta)
dy = dr2
dz = dr1 + dr3*np.sin(theta)

dx_test_if_dx_eqal_to_this = 2*np.pi/(N3*dq3)

# define FOV coodinates in xyz system. (r1 r2 r3 are meshgrids! )
x = r3*np.cos(theta)
y = r2
z = r1 + r3*np.sin(theta)

# define FOV of xyz (intensity values, not coordinates) ((ett intensitetsvärde för varje coordinat som är sparad i x y z meshgrids))
FOV_yzx = np.zeros((N2, N1, N3))
                    # y,z,x?

# size of FOV in xyz
FOV_size_z = dz*N1
FOV_size_y = dy*N2     # är dessa rätt?
FOV_size_x = dx*N3

#==============================================================================
# Create an object in the xyz system (a simple test object. to create real you have to (?) create a orthoganal XYZ system and then skew it to our xyz )
# 4 hörn?
#==============================================================================
# size of real object in [m]
x_obj = 300E-9 
y_obj = 100E-9 # är dessa xyz rätt?
z_obj = 500E-9

# how many pixels (object plane pixels) that the object sizes fits in 
Nx_obj = int(ceil( x_obj / dx))
Ny_obj = int(ceil( y_obj / dy)) # är dessa xyz rätt, alltså kallas de rätt?
Nz_obj = int(ceil( z_obj / dz))
#
## create object of size N_obj 
#obj = np.ones((Ny_obj,Nz_obj,Nx_obj))
#testtt= np.ones(())

# insert object into xyz FOV (centred)
#FOV_xyz[0,1,0] = 1

 #  Sned   y ,z , x
#FOV_yzx[int((N2-Ny_obj)/2) : int((N2-Ny_obj)/2) + Ny_obj, int((N1-Nz_obj)/2):int((N1-Nz_obj)/2) + Nz_obj, int((N3-Nx_obj)/2):int((N3-Nx_obj)/2) + Nx_obj ] = obj


#==============================================================================
# "Move" that simple object you created in xyz to r1r2r3
#==============================================================================
# TODO:
    


#==============================================================================
# Här börjar jag göra ett annat objekt, av "kristall"-form
#==============================================================================
#dy_obj_real = 30E-9   #verkliga avståndet mellan 2 lattice points
#dx_obj_real = 20E-9   #
#dz_obj_real = 20E-9   #
#
##hur många pixlar är det mellan 2 lattice points
#dy_obj = ceil(dy_real/pixel_objectPlane)    #this should be defined as a real distance instead of this pixel thing. otherwise the nbr of lattice points will change relative to z
#dx_obj = ceil(dx_real/pixel_objectPlane)
#dz_obj = ceil(dz_real/pixel_objectPlane)
#
#dyrange = [0, Ny_obj]           #borde vara relativt FOV!!!??
#dxrange = [0, Nx_obj]
#dzrange = [0, Nz_obj]
#
##this is so that i can do scattering plots
#lattice_points_z = np.ceil((dzrange[1] - dzrange[0]) / dz_obj)
#lattice_points_y = np.ceil((dyrange[1] - dyrange[0]) / dy_obj)
#lattice_points_x = np.ceil((dxrange[1] - dxrange[0]) / dx_obj)             
##total nbr of lattice points
#lattice_points = int(lattice_points_z * lattice_points_y *lattice_points_x)
#
#ycoor = np.zeros((lattice_points))
#xcoor = np.zeros((lattice_points))
#zcoor = np.zeros((lattice_points))
#
#index = 0 
#for layer in range(dzrange[0], dzrange[1], dz):
#    for row in range(dyrange[0], dyrange[1], dy):  
#        for col in range(dxrange[0], dxrange[1], dx):
#            zcoor[index] = layer
#            ycoor[index] = row    
#            xcoor[index] = col
#            obj[layer,row,col] = 1
#            index = index+1 
#
## insert the object into the FOV (centered in a 3 dim)
##FOV_l[int((Nz_det-Nz_obj)/2) : int((Nz_det-Nz_obj)/2) + Nz_obj, int((Nxy_det-Ny_obj)/2) : int((Nxy_det-Ny_obj)/2) + Ny_obj, int((Nxy_det-Nx_obj)/2) : int((Nxy_det-Nx_obj)/2) + Nx_obj ] = obj
#      
## FFT2 of the space in FOV at that distance z      
## this is what the pattern from the obj would look lika at that distance
##fft_FOV_l = fft.fftshift(fft.fft2(FOV_l))      
#
#transMatrix = np.zeros((4,4))
#theta1=theta
## skew FOV_l
#t =np.array([ [np.cos(theta1), 0,-np.sin(theta1),0],
#               [0,1,0, 0],[np.sin(theta1),0,np.cos(theta1),0],[0,0,0,1]])

def plot_3Dscatter():
#    #scatter 3 ritar ut en ring för varje punkt specificerad av vektorerna (x,y,z)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(r1,r2,r3, alpha=0.2, c='red', marker ='.')
    
#    ax.set_xlim3d(50,150)
#    ax.set_ylim3d(50,150)
#    ax.set_zlim3d(50,150)
#plot_3Dscatter()     

def plot_3Dscatter_R():
#    #scatter 3 ritar ut en ring för varje punkt specificerad av vektorerna (x,y,z)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(r3*1E6,r2*1E6, r1*1E6, alpha=0.2, c='red', marker ='+')
    #ax.plot(r1[:,:,0], r2[:,:,0], 'r+', zdir='r3', zs=1.5)
#    ax.set_xlim3d(50,150)
#    ax.set_ylim3d(50,150)
#    ax.set_zlim3d(50,150)
    plt.xlabel(' r1 [um]')
    plt.ylabel(' r2 [um]')
    plt.title('r1 r2 r3')
#plot_3Dscatter_R() 

def plot_3Dscatter_XYZ():
#    #scatter 3 ritar ut en ring för varje punkt specificerad av vektorerna (x,y,z)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x*1E6, y*1E6, z*1E6, alpha=0.2, c='blue', marker ='+')
    #ax.plot(r1[:,:,0], r2[:,:,0], 'r+', zdir='r3', zs=1.5)
#    ax.set_xlim3d(50,150)
#    ax.set_ylim3d(50,150)
#    ax.set_zlim3d(50,150)
    plt.xlabel(' x [um]')
    plt.ylabel(' y [um]')
    plt.title('x y z')
#plot_3Dscatter_XYZ() 

def plot_xyzFOV_object():
#    #scatter 3 ritar ut en ring för varje punkt specificerad av vektorerna (x,y,z)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x*1E6,y*1E6,z*1E6, c= FOV_yzx, marker ='o')
  
    #ax.plot(r1[:,:,0], r2[:,:,0], 'r+', zdir='r3', zs=1.5)
#    ax.set_xlim3d(50,150)
#    ax.set_ylim3d(50,150)
#    ax.set_zlim3d(50,150)
    plt.xlabel(' x [um]')
    plt.ylabel(' y [um]')
    ax.set_zlabel('z [um]')
    #plt.zlabel(' z [um]')
    plt.title('Simple object in xyz FOV.')
#plot_xyzFOV_object() 

#plt.figure()
#plt.subplot(221)
#pixel_objectPlane_nano = pixel_objectPlane *1E9
#plt.title('z=0, dx=%.2f nm '%pixel_objectPlane_nano)
#plt.imshow(FOV_l, extent=[0, FOV_real_size*1E6, 0, FOV_real_size*1E6])#, cmap='gray')  FOV_real_size
#plt.xlabel(' x [um]')
#plt.ylabel(' y [um]')
#
#plt.subplot(222)
#plt.title('z=%.2f m dx = 172 um'%z)
#plt.imshow(abs(fft_FOV_l), extent=[0, Nxy_det * pixel_det*1E3, 0, Nxy_det * pixel_det*1E3])#, cmap='gray')  FOV_real_size
#plt.xlabel(' x [mm]')
#plt.ylabel(' y [mm]')
#
#plt.savefig('savefig\diff_sim_crystal_%d'%z, bbox_inches='tight')
