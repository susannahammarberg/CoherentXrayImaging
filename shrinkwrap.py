# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:04:38 2017

@author: HonkyT
"""
#from numpy import fft # om man skriver såhär begöver man itne använda np.
# men eftersom jag behöver fler saker från nupy och jag redan skrivit np på alla gör jag inte det.
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from create2Dgaussian import create2Dgaussian  #import all functions
from scipy import misc
from math import ceil

data = np.genfromtxt("diffractionPattern.txt", dtype='float')
#image = misc.imread('P.png',flatten=True)
#p = misc.imread('star.bmp',flatten=True)
#absF = abs(np.fft.fftshift(np.fft.fft2(p)))
########### Test att använda fft av kristall som inputi shrinkwrap
crystal = np.zeros((201,201), dtype= np.int32)

dx=1
for row in range(60,140,dx):
    for col in range(60,140,dx):
        crystal[row,col] = 1

crystal_fourier = np.fft.fftshift(np.fft.fft2(crystal))
absF = abs(crystal_fourier)

#############Test att använda bild som är skjuvad coord
#image = misc.imread('star.bmp',flatten=True)
##image = misc.imread('P.png',flatten=True)
#data = abs(np.fft.fftshift(np.fft.fft2(image)))
#
#theta = 0*np.pi / 180 #rad
#r3 = 1 + 1/np.cos(theta)    #antal pixlar som motsvarar 1 pixel i xled i xz systemet
#r1 = 1 + 1/ np.sin(theta)   
#
#image_skewed = np.zeros((image.shape[0], ceil((image.shape[1]/np.cos(theta)))+ 50  ))
#
#for i in range(0,image.shape[0]):
#    for j in range(0,image.shape[1]):
#        xs = ceil(j / (1 + np.tan(theta)**2) + i*np.tan(theta)/ ( 1 + np.tan(theta)**2)  )
#        #np.disp(xs)
#        ys = i
#        image_skewed[ys,xs] = image[i,j] 
#        #image_skewed = 
#
## cut ot hte skewed image so that it is ccentered
#image_skewed = image_skewed[:,0:120] # Hamnar rätt om man klipper den rätt
#image_skewed = image_skewed/image_skewed.max()
#fft_image_skewed = abs(np.fft.fftshift(np.fft.fft2(image_skewed)))
#absF = fft_image_skewed 

#fattar inte hur man skriver det som function på ett smidigt sätt.
# går ju inte att kolla på variablerna om man kör som funktion
#def Shrinkwrap(diffPattern):
    # the intensity of the diffraction pattern is F^2
    
#absF = pow(data,0.5)

# parameters
beta = 0.9               # HIO feedback parameter
thresholdMask1 = 0.04;   # Initial threshold for mask
thresholdMask2 = 0.2;    # threshold for mask 
sigma = 3;               # initial value of gaussian width
c = absF.shape[0]/2      # center of Gaussian

n = 100                   #number of iterations

# initialize vector for error calculation
err2 = np.zeros(shape=(1,n))

# g represents the unknown object. Create initial guess for g using 
# the autocorrelation function of the object. Start by Fourier 
# transform the diffraction pattern:
g = np.fft.fftshift(np.fft.ifft2(absF)) #finns några val som symmetric?
# create autocorrelation function of g
#kanske finns än funtion i python som gör det här
g = abs(np.fft.fftshift(np.fft.ifft2(np.fft.fft2(g)*np.conj(np.fft.fft2(g)))))/(g.shape[0]*g.shape[1])
# create a logical matrix that is the Shrinkwrap support
support = g > thresholdMask1*g.max()
# define iteration counter
k = 0
# Start of iterations 
while k<n: 
    # every 20th iteration, update the support
    if k%20==0 and k>0:   # VARFÖR SKULLE K VA NEG?
        # call function create2Dgaussian
        gauss2D = create2Dgaussian(sigma,c,absF.shape[0],absF.shape[1])
        # calculate the convolution the absolute value of the object wave 
        # with a gaussian
        # finns en djungel av olika laternativ hur man convolverar
        support = signal.convolve2d(gauss2D,abs(g),'same')
        # Create a logical matrix that is the Shrinkwrap support
        support = support > thresholdMask2*support.max()                        #SÄKERT ATT DEN HITTAR 2D MAXET?      
        # reduce sigma with 1% every 20th iteration until a limit of 1.5
        if sigma >= 1.5:
            sigma = sigma*0.99
        
    # STEP 1: Fourier transform of g(x),  G(u) :
    G = np.fft.fftshift(np.fft.fft2(g))  
    
    # STEP 2: make |G| confirm with |F|
    Gprime = absF*np.exp(1j*np.angle(G))
    
    # STEP 3: inverse Fourier transform      
    gprime = np.fft.ifft2(np.fft.ifftshift(Gprime))
    
    # STEP 4: See to that gprime satisfies its constraints:   
            
    # create inverse of the support (compuational aid)    
    support_inv = np.logical_not(support)
   
    # update g(l,j) outside the support
    gprime = gprime*support + g*support_inv - beta*gprime*support_inv
    
    # update g'(l,j) inside the support in one of 
    # the two following ways:
    # gprime(l,j) = g(l,j)(in matrix multiplication i.e.:)
    #gprime = gprime*support_inv + g*support
    #or    
    # gprime[l,j] = gprime[l,j]    (i.e. not at all)
        
        
        
#    for l in range(0,absF.shape[0]):
#       for j in range(0,absF.shape[1]):         ## ska DET VARA RANGE 0 ÄR FÖRSTA INDEXEDT 0?
#           if support[l,j] == 0:   
#               # update g(l,j) outside the support
#               gprime[l,j] = g[l,j] - beta * gprime[l,j]
#               
#           else:
#               # update g'(l,j) inside the support in one of 
#               # the two following ways:
#               #gprime(l,j) = g(l,j);
#               gprime[l,j] = gprime[l,j]

    # overwrite g with result from iteration
    g = gprime
    # set all negative values of g to zero
    g[g<0] = 0  

    # för att lätt kunna titta på g
    realg = g.real
    
    err2[0][k] = sum(sum( (abs(abs(np.fft.fftshift(np.fft.fft2(((g*np.conj(g))))))) - absF)**2)) / sum(sum(absF**2))

    k = k+1
    
    
plt.figure()    
plt.imshow(abs(g), cmap='gray')
plt.title('reconstructed object')

plt.figure()
plt.imshow(crystal, cmap='gray')
plt.title('original object')