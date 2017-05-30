# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 11:19:26 2017

@author: Susanna Hammarberg + snott från JWallentin
Simulation based on the orientasion in Berenguer et al.

"""
from mpl_toolkits.mplot3d import Axes3D
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt
from math import ceil
from numpy import fft

# Välj dessa parameterar så att pixelstorleken blir så lätt som möjligt
# så att de liksom hamnar på griden. om det går. 
# Simulation parameters
L = 2.8              # Object-detector distance, m; 0.4357 ; 0.5849 ; 0.8919 ; 1.4804
# pixel size on detector
pixel_det = 22E-6   # Pixel ctr to ctr distance (w) [m]
energy = 9.5        # keV
wavelength = 1.23984E-9/energy
kLength = 2*np.pi/wavelength      # 1/m   needed?
theta = 10*np.pi / 180
dtheta = 0.001*np.pi/180
# pixel size real space 
#deltaS = wavelength * L / ( N * pixel_det)

# create matrix to hold diffraction patterns
#diffSet=np.zeros((nbr_scans, 195, 487))

#def pixel_sizes():
        # number of pixels in detector space (the once that are used)
Nq1 = 35
Nq2 = 29
Nq3 = 11 #number of steps on rocking curve. q3=0 should be in the middle of q3

dq1 = 1E-9*np.pi * pixel_det / ( wavelength * L) 
dq2 = dq1

# length of scattering vector (hur får man det?)length of k =2pi/lambda
q_abs = 26.7E9        #1/m
dq3 = 1E-9 * q_abs*dtheta

## Create r1 r2 r3 space
## according to Jespers pdf
dr1 = 2*np.pi / Nq1 * dq1 * np.cos(theta)
dr2 = 2*np.pi / Nq2 * dq2 * np.cos(theta)
dr3 = 2*np.pi / Nq3 * dq3 * np.cos(theta)

# create realspace 
# pixel sizes
dx = dr3*np.cos(theta)
dy = dr2
dz = dr1 + dr3*np.sin(theta)

#
#plt.close("all")
#def create_crystal_object():
#    # Create simulation object with phase
#    #####################################
#    crystal3D = np.zeros((201,201,201), dtype= np.int32)
#    crystal3D_phase = np.zeros((201,201,201))
#    dx=3
#    for row in range(60,140,dx):
#        for col in range(80,120,dx):
#            for time in range(90,110,dx):
#                crystal3D[row,col,time] = 1
#                # add phase
#                if np.random.randint(2) == 1:
#                    crystal3D_phase[row,col,time] = -np.pi/1
#                else:
#                    crystal3D_phase[row,col,time] = np.pi/1
#     
#    # crystal =    
#    return dx
#                
## create meshgrid for r1,r2,r3
## number of pixels in detector space (the once that are used)
#Nq1 = 35
#Nq2 = 29
#Nq3 = 11 #number of steps on rocking curve. q3=0 should be in the middle of q3
#
## inte samma som J MEN detta är ju inte heller rätt!
## distances in detector plane to object plane
#dq1 = L * wavelength /( Nq1 * pixel_det) 
#dq2 = L * wavelength /( Nq2 * pixel_det)
#
## length of scattering vector (hur får man det?)length of k =2pi/lambda
#q_abs = 26.7E9 
#dq3 = q_abs*dtheta      *1E-12
#
## inte samma som J
#q1_linspace = np.linspace( - dq1*Nq1 /2, dq1*Nq1 /2, Nq1 )
#q2_linspace = np.linspace( - dq2*Nq2 /2, dq2*Nq2 /2, Nq2 )
#q3_linspace = np.linspace( - dq3*Nq3 /2, dq3*Nq3 /2, Nq3 )
#
#q1, q2, q3 = np.meshgrid( q1_linspace, q2_linspace, q3_linspace )
#
## Create r1 r2 r3 space
## according to Jespers pdf
#dr1 = 2*np.pi / Nq1 * dq1 * np.cos(theta)
#dr2 = 2*np.pi / Nq2 * dq2 * np.cos(theta)
#dr3 = 2*np.pi / Nq3 * dq3 * np.cos(theta)
#
## J uses floor . why?
#r1_linspace = np.linspace( - dr1*Nq1 /2, dr1*Nq1 /2, Nq1 )
#r2_linspace = np.linspace( - dr2*Nq2 /2, dr2*Nq2 /2, Nq2 )
#r3_linspace = np.linspace( - dr3*Nq3 /2, dr3*Nq3 /2, Nq3 )
#
#r1, r2, r3 = np.meshgrid(r1_linspace, r2_linspace, r3_linspace)
#
## create realspace (!? ) dessa tre 
## pixel sizes
#dx = dr3*np.cos(theta)
#dy = dr2
#dz = dr1 + dr3*np.sin(theta)
#
## x y z FOV: (nöjer mig med detta sålänge)
## dessa tre vektorer måste ju vara ortogonala
#x = r3*np.cos(theta)
#y = r2
#z = r1 + r3*np.sin(theta)
#
#
## test nonsence
#t1 = [1, 2, 3 ]
#t2 = [1 ,3 ,9 ]
#t3 = [1 ,1, 1 ]
#def plot_crystal3D():
#    #scatter 3 ritar ut en ring för varje punkt specificerad av vektorerna (x,y,z)
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    ax.scatter(z,y,x, c='red')#, marker ='.'
#    
#    
#    plt.xlabel(' x')
#    plt.ylabel(' y')
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    ax.scatter(r1,r2,r3, c='red', marker ='.')
#
##    plt.subplot(222)
##    plt.imshow(abs(crystal3D[102,:,:]), cmap='gray')
##    #plt.title('xz cut')
##    plt.xlabel(' z')
##    plt.ylabel(' x')
##    plt.colorbar()
##    
##    plt.subplot(223)        #x
##    plt.imshow(abs(crystal3D[:,101,:]), cmap='gray')
##    #plt.title('yz cut')
##    plt.xlabel(' z')
##    plt.ylabel(' y')
##    plt.colorbar()
#    
#    
#plot_crystal3D()    
#
## create object in x,y,z
#


#x = np.linspace(0,256,257)
#y = np.linspace(0,256,257)
#z = np.linspace(0,256,257)

#r1 =  x/np.cos(theta)
#r2 = y
#r3 = z - x/np.tan(theta)
#
## kolla hur J gör meshgrid till ett snet coordinatsystem
#
##for row in range(())
#x, y  = np.meshgrid(x,y)    # gör 2 matriser med coordinater

