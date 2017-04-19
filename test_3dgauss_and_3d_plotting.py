# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 15:50:57 2017

@author: HonkyT
"""
import numpy as np
from numpy import fft
from scipy import signal
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from create3Dgaussian import create3Dgaussian
from math import ceil
from scipy import misc

plt.close("all")


# create and plot gaussian
############################
def gaussian():
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

# Create crystal: skewed and unskewed
# create FFTs of these. plot
################################################
image = misc.imread('star.bmp',flatten=True)
#image = misc.imread('P.png',flatten=True)
data = abs(fft.fftshift(fft.fft2(image)))

theta = 10*np.pi / 180 #rad
r3 = 1 + 1/np.cos(theta)    #antal pixlar som motsvarar 1 pixel i xled i xz systemet
r1 = 1 + 1/ np.sin(theta)   

image_skewed = np.zeros((image.shape[0], ceil((image.shape[1]/np.cos(theta)))+ 50  ))

for i in range(0,image.shape[0]):
    for j in range(0,image.shape[1]):
        xs = ceil(j / (1 + np.tan(theta)**2) + i*np.tan(theta)/ ( 1 + np.tan(theta)**2)  )
        #np.disp(xs)
        ys = i
        image_skewed[ys,xs] = image[i,j] 
        #image_skewed = 

fft_image = fft.fftshift(fft.fft2(image))        
fft_image_skewed = fft.fftshift(fft.fft2(image_skewed))

# look at crystal and their ffts
crystal = np.zeros((201,201), dtype= np.int32)
crystal3D = np.zeros((201,201,201), dtype= np.int32)
crystal3D_filter = np.zeros((201,201,201), dtype= np.int32)
crystal_filter2D=np.zeros((201,201), dtype= np.int32)
dx=3
dx_filter = 1

for row in range(60,140,dx):
    for col in range(60,140,dx):
        crystal[row,col] = 1

for row in range(60,140,dx):
    for col in range(80,120,dx):
        for time in range(90,110,dx):
            crystal3D[row,col,time] = 1
           
            
# construct filters            
for row in range(92,108,dx_filter):
    for col in range(92,108,dx_filter):
        crystal_filter2D[row,col] = 1
        
for row in range(80,120,dx_filter):
    for col in range(80,120,dx_filter):
        for time in range(80,120,dx_filter):
            crystal3D_filter[row,col,time] = 1
            
crystal_skewed = np.zeros((crystal.shape[0], ceil((crystal.shape[1]/np.cos(theta)))+ 50  ))

for i in range(0,crystal.shape[0]):
    for j in range(0,crystal.shape[1]):
        xs = ceil(j / (1 + np.tan(theta)**2) + i*np.tan(theta)/ ( 1 + np.tan(theta)**2)  )
        #np.disp(xs)
        ys = i
        crystal_skewed[ys,xs] = crystal[i,j] 
        #image_skewed = 


       
fft_crystal_skewed = fft.fftshift(fft.fft2(crystal_skewed))

crystal3D_fourier = fft.fftshift(fft.fftn(crystal3D))
crystal_fourier = fft.fftshift(fft.fft2(crystal))

diffPattern3D = abs(crystal3D_fourier)**2


# keep only one bragg peak of the FFT3
filtered_crystal3D_fourier = crystal3D_filter * crystal3D_fourier
# ifft back to real space
filtered_real = abs(fft.ifftn(fft.ifftshift(filtered_crystal3D_fourier)))

#projectionMatrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0]])
#projectedCrystal = projectionMatrix*crystal3D

# plot  crystal
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

#plt.figure()
#plt.imshow(projectedCrystal)
#plt.title('Orthoganl projection of 3D crystal into x-y-plane')
#plt.xlabel(' x')
#plt.ylabel(' y')
#savefig('test.jpg')
########################plot 3D crystal in all plane cuts 
plt.figure()
plt.subplot(221)
plt.imshow(crystal3D[:,:,102], cmap='gray')
#plt.title('xyplane cut')
plt.xlabel(' x')
plt.ylabel(' y')


plt.subplot(222)
plt.imshow(crystal3D[102,:,:], cmap='gray')
#plt.title('xz cut')
plt.xlabel(' z')
plt.ylabel(' x')


plt.subplot(223)        #x
plt.imshow(crystal3D[:,101,:], cmap='gray')
#plt.title('yz cut')
plt.xlabel(' z')
plt.ylabel(' y')


# plot filtered 2D retrived crystal2D
plt.figure
fig,ax = plt.subplots(1)
plt.imshow(abs((fft.ifft2(fft.ifftshift(test)))))

########################plot 3D crystal in all plane cuts in reciprocal space
plt.figure()
plt.subplot(221)
plt.imshow(np.log10(abs(crystal3D_fourier[:,:,100])),  cmap='gray')
#plt.title('xy plane cut')
plt.xlabel(' x')
plt.ylabel(' y')
plt.colorbar()

plt.subplot(222)
plt.imshow(np.log10(abs(crystal3D_fourier[100,:,:])),  cmap='gray')
#plt.title('xz plane ')
plt.xlabel(' z')
plt.ylabel(' x')    #rätt
plt.colorbar()

plt.subplot(223)
plt.imshow(np.log10(abs(crystal3D_fourier[:,102,:])), cmap='gray')
#plt.title('yz plane')
plt.xlabel(' z')
plt.ylabel(' y')
plt.colorbar()

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
    
    #plt.axis.set_xticklabels([])
    
    plt.subplot(223)
    plt.imshow(np.log10(abs(filtered_crystal3D_fourier[:,102,:])), cmap='gray')
    #plt.title('yz plane')
    plt.xlabel(' z')
    plt.ylabel(' y')
    plt.colorbar()

plot_filtered_crystalFFt()

def plot_filtered_realspace():
    plt.figure()
    plt.subplot(221)
    plt.imshow(filtered_real[:,:,102],  cmap='gray')
    #plt.title('xy plane cut')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.colorbar()
    
    plt.subplot(222)
    plt.imshow(filtered_real[100,:,:],  cmap='gray')
    #plt.title('xz plane ')
    plt.xlabel(' z')
    plt.ylabel(' x')    #rätt
    plt.colorbar()
    
    #plt.axis.set_xticklabels([])
    
    plt.subplot(223)
    plt.imshow(filtered_real[:,102,:], cmap='gray')
    #plt.title('yz plane')
    plt.xlabel(' z')
    plt.ylabel(' y')
    plt.colorbar()

    
plot_filtered_realspace()
    
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
