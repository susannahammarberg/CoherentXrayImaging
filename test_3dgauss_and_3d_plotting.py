# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 15:50:57 2017

@author: Susanna Hammarberg
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

# load image. look at FFT2. Skew image. Look att FFT2.
def image():
    image = misc.imread('star.bmp',flatten=True)
    circle = misc.imread('circle.png',flatten=True)
    low_values_indices = circle < 200  # Where values are low
    circle[low_values_indices] = 0  # All low values set to 0
    high_values_indices = circle > 0
    circle[high_values_indices] = 1
    #image = misc.imread('P.png',flatten=True)
    #data = abs(fft.fftshift(fft.fft2(image)))
    
    theta = 10*np.pi / 180 #rad
    # nä detta stämmer väl inte
    r3 = 1 + 1/np.cos(theta)    #antal pixlar som motsvarar 1 pixel i xled i xz systemet
    r1 = 1 + 1/ np.sin(theta) 
    return 0

def skew_image():
    image_skewed = np.zeros((image.shape[0], ceil((image.shape[1]/np.cos(theta)))+ 50  ))
    
    for i in range(0,image.shape[0]):
        for j in range(0,image.shape[1]):
            xs = ceil(j / (1 + np.tan(theta)**2) + i*np.tan(theta)/ ( 1 + np.tan(theta)**2)  )
            #np.disp(xs)
            ys = i
            image_skewed[ys,xs] = image[i,j] 
            
    fft_image_skewed = fft.fftshift(fft.fft2(image_skewed))
    return image_skewed, fft_image_skewed
# image in reciprocal space
#fft_image = fft.fftshift(fft.fft2(image))


# Create crystal: skewed and unskewed. 2D and 3D.
# create FFTs of these. Filter out one peak.
################################################
Nx = 350
Ny = 350
Nz = 350
crystal = np.zeros((201,201), dtype= np.int32)
crystal3D = np.zeros((Nz,Ny,Nx))
crystal3D_filter = np.zeros((Nz,Ny,Nx), dtype= np.int32)
# coordinates of crystal points (for 3D scattering plot)
ycoor = np.zeros((6000))
xcoor = np.zeros((6000))
zcoor = np.zeros((6000))#2646
crystal3D_phase = np.zeros((Nz,Ny,Nx))

crystal_filter2D=np.zeros((201,201), dtype= np.int32)
dx = 3
dx_filter = 1

for row in range(60,140,dx):
    for col in range(60,140,dx):
        crystal[row,col] = 1

index=0
#==============================================================================
# for row in range(20,180,dx):  
#     for col in range(80,120,dx):
#         for time in range(90,110,dx): 
#==============================================================================
for row in range(80,120,dx):  
    for col in range(80,120,dx):
        for time in range(90,110,dx):           
            ycoor[index] = row    #(these are for 3Dscattering plot)
            xcoor[index] = col
            zcoor[index] = time
            
            crystal3D[row,col,time] = 1
            # add phase
            if np.random.randint(2) == 1:
                crystal3D_phase[row,col,time] = -np.pi/1
            else:
                crystal3D_phase[row,col,time] = np.pi/1
            index = index+1 


# make an object of both phase and amplitude            
crystal3D = crystal3D * np.exp(1j*crystal3D_phase)

def pad_3D_object(obj, Nx, Ny, Nz): 
    padded_object = np.zeros((Nz, Ny, Nx), dtype=np.complex64)
    x = (Nx - obj.shape[2]) / 2
    y = (Ny - obj.shape[1]) / 2 
    z = (Ny - obj.shape[0]) / 2 
    padded_object[z: z+ obj.shape[0] , y: y + obj.shape[1], x: x+ obj.shape[2]] = obj
    
    return padded_object      


# warning: if you use the padding the 3D scatter plot will not show you that
# innessessary or easier than to redefine the matris from the start?
#crystal3D = pad_3D_object(crystal3D, 350, 350, 350)#   350

                           

# construct 2D and 3D filters            
for row in range(92,108,dx_filter):
    for col in range(92,108,dx_filter):
        crystal_filter2D[row,col] = 1
        
for row in range(80,120,dx_filter):
    for col in range(80,120,dx_filter):
        for time in range(80,120,dx_filter):
            crystal3D_filter[row,col,time] = 1
    


def skewed_crystal_diffPatterns(crystal, theta):    
    dx_filter = 1            
    crystal_skewed = np.zeros((crystal.shape[0], np.floor((crystal.shape[1]/np.cos(theta))) + 50  ))
    crystal_filter2D_skewed =  np.zeros(( crystal_skewed.shape[0], crystal_skewed.shape[1] ))
    # simulate skewed diffraction patterns
    theta = 10*np.pi/180
    
    for i in range(0,crystal.shape[0]-6):
        for j in range(0,crystal.shape[1]-6):
            xs = ceil(j / (1 + np.tan(theta)**2) + i*np.tan(theta)/ ( 1 + np.tan(theta)**2)  )
            ys = i
            crystal_skewed[ys,xs] = crystal[i,j] 
    
    for row in range(80,150,dx_filter):
        for col in range(80,160,dx_filter):            
            crystal_filter2D_skewed[row,col] = 1
 
    diffPattern = abs( fft.fftshift(fft.fft2(crystal_skewed)))**2
    return (diffPattern, crystal_filter2D_skewed)

#diffPattern_skewed, crystal_filter2D_skewed = skewed_crystal_diffPatterns(crystal, theta )


# call shrinkwrap for filtered skewed crystal
# yel = shrinkwrap(diffPattern_skewed*crystal_filter2D_skewed)



sum_y_crystal3D = np.sum(crystal3D,axis=0)
#sum_y_crystal3D_fourier = fft.fftshift(fft.fft2(sum_y_crystal3D))
#sumtest = abs(sum_y_crystal3D_fourier)**2
#sum_y_crystal3D = shrinkwrap(sumtest)

sum_y_crystal3D_angle = np.zeros((201,201,201))

for row in range(0,201):
    for col in range(0,201):
        sum_y_crystal3D_angle[row,col] = sum_y_crystal3D[row,col] * np.cos(theta)

test_diffpattern = abs(fft.fftshift(fft.fft2(sum_y_crystal3D_angle)))**2



crystal3D_fourier = fft.fftshift(fft.fftn(crystal3D))
crystal_fourier = fft.fftshift(fft.fft2(crystal))


# keep only one bragg peak of the FFT3
#filtered_crystal3D_fourier = crystal3D_filter * crystal3D_fourier

#diffPattern3D = abs(filtered_crystal3D_fourier)**2

# ifft back to real space
#filtered_real = fft.ifftn(fft.ifftshift(filtered_crystal3D_fourier))

# take aboluter value and do ifft back to real space
#filtered_real = ( fft.ifftn(fft.ifftshift(abs(filtered_crystal3D_fourier))))

# call Shrinkwrap
#filtered_real = shrinkwrap3D(diffPattern3D)


#diffPattern2Dy = np.sum(diffPattern3D, axis=0)
#diffPattern2Dx = np.sum(diffPattern3D, axis=1)
#diffPattern2Dz = np.sum(diffPattern3D, axis=2)


# call shrinkwrap 2D
#retrieved_y = shrinkwrap(diffPattern2Dy)
#retrieved_x = shrinkwrap(diffPattern2Dx)

#########################################
# Plot functions
#########################################

def plot_3Dscatter():
#    #scatter 3 ritar ut en ring för varje punkt specificerad av vektorerna (x,y,z)
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    ax.scatter(z,y,x, c='red')#, marker ='.'
#    plt.xlabel(' x')
#    plt.ylabel(' y')
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xcoor,ycoor,zcoor, alpha=0.2, c='red', marker ='.')
    ax.set_xlim3d(50,150)
    ax.set_ylim3d(50,150)
    ax.set_zlim3d(50,150)
plot_3Dscatter()

def plot_skewed_2Dcrystal():
    plt.figure()
    plt.imshow(diffPattern_skewed, cmap='gray')
    plt.colorbar()
    
    plt.figure()
    plt.imshow(crystal_filter2D_skewed, cmap='gray')
    plt.colorbar()
    
    plt.figure()
    plt.imshow(diffPattern_skewed*crystal_filter2D_skewed, cmap='gray')
    plt.colorbar()
#plot_skewed_2Dcrystal()

def plot2dsummed(retrieved_x):
    plt.figure()
    plt.subplot(211)
    plt.imshow(abs(retrieved_x), cmap='gray')
    plt.suptitle('Crystal amplitude (above) and phase (below) (xdim) ')
    plt.colorbar()
    plt.subplot(212)
    plt.imshow(np.angle(retrieved_x), cmap='gray')
    plt.colorbar()

#plot2dsummed()
#plot2dsummed(yel)

def plot_crystal2D():
        
    plt.figure()
    plt.subplot(221)
    plt.imshow(crystal, cmap='gray')
    plt.title('crystal')
    plt.axis('off')
    
    plt.subplot(223)
    test=crystal_filter2D*crystal_fourier
    plt.imshow(np.log10(abs(test)), cmap='gray')
    plt.title('log10 |FFT2| of crystal')
    plt.axis('off')
    
    plt.subplot(222)
    plt.imshow(crystal_skewed, cmap='gray')
    plt.xlabel(' x_s')
    plt.ylabel(' y_s')
    plt.title('crystal in skewed coordinates')
    plt.axis('off')
       
    plt.subplot(224)
    plt.imshow(np.log10(abs(fft_crystal_skewed)), cmap='gray')
    plt.xlabel(' 1/x_s')
    plt.ylabel(' 1/y_s')
    plt.title('log10|FFT| of crystal in skewed coordinates')
    plt.axis('off')

    # plot filtered 2D retrived crystal2D
    plt.figure
    fig,ax = plt.subplots(1)
    plt.imshow(abs((fft.ifft2(fft.ifftshift(test)))))

#plot_crystal2D()


#plt.figure()
#plt.imshow(projectedCrystal)
#plt.title('Orthoganl projection of 3D crystal into x-y-plane')
#plt.xlabel(' x')
#plt.ylabel(' y')
#savefig('test.jpg')
########################plot 3D crystal in all plane cuts 
def plot_crystal3D():
    plt.figure()
    plt.subplot(221)
    plt.imshow(abs(crystal3D[:,:,176]), cmap='gray')
    #plt.title('xyplane cut')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.subplot(222)
    plt.imshow(abs(crystal3D[175,:,:]), cmap='gray')
    #plt.title('xz cut')
    plt.xlabel(' z')
    plt.ylabel(' x')
    plt.colorbar()
    
    plt.subplot(223)        #x
    plt.imshow(abs(crystal3D[:,175,:]), cmap='gray')
    #plt.title('yz cut')
    plt.xlabel(' z')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.suptitle('Plane cuts of crystal amplitude')
        
    plt.figure()   #Phase
    plt.subplot(221)
    plt.imshow(np.angle(crystal3D[:,:,176]), cmap='gray')
    #plt.title('xyplane cut')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.subplot(222)
    plt.imshow(np.angle(crystal3D[175,:,:]), cmap='gray')
    #plt.title('xz cut')
    plt.xlabel(' z')
    plt.ylabel(' x')
    plt.colorbar()
    
    plt.subplot(223)        #x
    plt.imshow(np.angle(crystal3D[:,175,:]), cmap='gray')
    #plt.title('yz cut')
    plt.xlabel(' z')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.suptitle('Plane cuts of crystal phase')
    
plot_crystal3D()    


########################plot 3D crystal in all plane cuts in reciprocal space
def plot_crystal3D_reciprocal():
    plt.figure()
    #plt.title('hej')
    plt.subplot(221)
    plt.imshow(np.log10(abs(crystal3D_fourier[:,:,174])),  cmap='gray')
    #plt.title('xy plane cut')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.subplot(222)
    plt.imshow(np.log10(abs(crystal3D_fourier[173,:,:])),  cmap='gray')
    #plt.title('xz plane ')
    plt.xlabel(' z')
    plt.ylabel(' x')    #rätt
    plt.colorbar()
    
    plt.subplot(223)
    plt.imshow(np.log10(abs(crystal3D_fourier[:,173,:])), cmap='gray')
    #plt.title('yz plane')
    plt.xlabel(' z')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.suptitle('Plane cuts of crystal in reciprocal space')

plot_crystal3D_reciprocal()

def plot_filtered_crystalFFt():
    plt.figure()
    plt.subplot(221)
    plt.imshow(np.log10(abs(filtered_crystal3D_fourier[:,:,100])),  cmap='gray')
    #plt.title('xy plane cut')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.subplot(222)
    plt.imshow(np.log10(abs(filtered_crystal3D_fourier[100,:,:])),  cmap='gray')
    #plt.title('xz plane ')
    plt.xlabel(' z')
    plt.ylabel(' x')    #rätt
    plt.colorbar()

    plt.subplot(223)
    plt.imshow(np.log10(abs(filtered_crystal3D_fourier[:,102,:])), cmap='gray')
    #plt.title('yz plane')
    plt.xlabel(' z')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.suptitle('Plane cuts of filtered crystal in reciprocal space')

#plot_filtered_crystalFFt()

def plot_filtered_realspace():
    
    plt.figure()
    plt.subplot(221)
    plt.imshow(abs(filtered_real[:,:,100]),  cmap='gray')
    #plt.title('xy plane cut')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.subplot(222)
    plt.imshow(abs(filtered_real[100,:,:]),  cmap='gray')
    #plt.title('xz plane ')
    plt.xlabel(' z')
    plt.ylabel(' x')    #rätt
    plt.colorbar()
    
    plt.subplot(223)
    plt.imshow(abs(filtered_real[:,102,:]), cmap='gray')
    #plt.title('yz plane')
    plt.xlabel(' z')
    plt.ylabel(' y')
    plt.colorbar()
    plt.suptitle('Plane cuts of \'retrieved\' crystal amplitude')

    
    plt.figure()
    plt.subplot(221)
    plt.imshow(np.angle(filtered_real[:,:,100]),  cmap='gray')
    #plt.title('xy plane cut')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.subplot(222)
    plt.imshow((np.angle(filtered_real[100,:,:])),  cmap='gray')
    #plt.title('xz plane ')
    plt.xlabel(' z')
    plt.ylabel(' x')    #rätt
    plt.colorbar()
    
    plt.subplot(223)
    plt.imshow(np.angle(filtered_real[:,102,:]), cmap='gray')
    #plt.title('yz plane')
    plt.xlabel(' z')
    plt.ylabel(' y')
    plt.colorbar()
    plt.suptitle('Plane cuts of \'retrieved\' crystal phase')
    
#plot_filtered_realspace()
    
# plot image and skewed image
def plot_star():
    
    plt.figure()
    plt.imshow(image)
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.title('Image in normal coordinates')
    
    plt.figure()
    plt.imshow(image_skewed)
    plt.xlabel(' x_s')
    plt.ylabel(' y_s')
    plt.title('Image in skewed coordinates')
    
    
    plt.figure()
    plt.imshow(np.log10(abs(fft_image)))
    plt.xlabel(' 1/x')
    plt.ylabel(' 1/y')
    plt.title('log10 abs of FFT of image')
    
    plt.figure()
    plt.imshow(np.log10(abs(fft_image_skewed)))
    plt.xlabel(' 1/x_s')
    plt.ylabel(' 1/y_s')
    plt.title('log10 abs of FFT of image in skewed coordinates')
    
    return 0


# create and plot gaussian
############################
def create_plot_gaussian():
    absF= np.ones((20,20,20))
    sigma = 3;               # initial value of gaussian width
    
    c = absF.shape[0]/2 #varför blir det en float?      # center of Gaussian
    
    gauss3D = create3Dgaussian(sigma,c,absF.shape[0],absF.shape[1],absF.shape[2])

    plt.figure()
    plt.subplot(221)
    plt.imshow(gauss3D[:,:,5], vmin=0, vmax=1,  cmap='gray')
    plt.title('3D Gauss[:,:,5]')
    plt.axis('off')
    
    plt.subplot(222)
    plt.imshow(gauss3D[:,:,10], vmin=0, vmax=1, cmap='gray')
    plt.title('3D Gauss[:,:,10]')
    plt.axis('off')
    
    plt.subplot(223)
    plt.imshow(gauss3D[1,:,:], vmin=0, vmax=1,  cmap='gray')
    plt.title('3D Gauss[1,:,:]')
    plt.axis('off')
    
    plt.subplot(224)
    plt.imshow(gauss3D[10,:,:], vmin=0, vmax=1,  cmap='gray')
    plt.title('3D Gauss[10,:,:]')
    plt.axis('off')
    return 0

#gaussian()

#import matplotlib.transforms as mtransforms
#
#def get_image():
#    from scipy import misc
#    Z = misc.imread('P.jpg')
#    return Z
#
## Get image
#fig, ax = plt.subplots(1,1)
#Z = get_image()
#
## image skew
#im = ax.imshow(Z, interpolation='none', origin='lower',
#                 extent=[-2, 4, -3, 2], clip_on=True)
#im._image_skew_coordinate = (3, -2)
#
#plt.show()
#

#TODO  create a polyhedron  
#def create_polyhedron():        
#    #dpol = 3            
#    #mini = 96
#    #poly = np.zeros((201,201,201))
#    #for row in range(mini+4,104,dpol):
#    #    for col in range(mini,120,dpol):
#    #        for time in range(mini,110,dpol):
#    #            poly[row,col,time] = 1
#    #            # add phase
#    #   # mini = mini-1        
#    return 0

    

