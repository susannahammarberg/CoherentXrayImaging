# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:04:38 2017

@author: HonkyT
"""
#from IPython import get_ipython
#get_ipython().magic('reset -sf')
#

#from numpy import fft # om man skriver såhär begöver man itne använda np.
# men eftersom jag behöver fler saker från nupy och jag redan skrivit np på alla gör jag inte det.
import numpy as np
from numpy import fft
from scipy import signal
from scipy.signal import fftconvolve
import matplotlib.pyplot as plt
from create3Dgaussian import create3Dgaussian  #import all functions
from scipy import misc
from math import ceil

def crystal():
    # look at crystal and their ffts
    crystal3D = np.zeros((201,201,201), dtype= np.int32)
    crystal3D_fourier = np.zeros((201,201,201), dtype= np.complex64)
    
    dx = 1
    
    for row in range(60,140,dx):
        for col in range(80,120,dx):
            for time in range(90,110,dx):
                crystal3D[row,col,time] = 1
    
    crystal3D_fourier = fft.fftshift(fft.fftn(crystal3D))
    #del crystal3D
    diffPattern3D = (abs(crystal3D_fourier)**2)
    del crystal3D_fourier
    return diffPattern3D

diffPattern3D = crystal()

def shrinkwrap3D(diffPattern3D):
      
    absF = np.zeros((diffPattern3D.shape[0], diffPattern3D.shape[1], diffPattern3D.shape[2]),dtype=np.complex64)  
    #gprime = np.zeros((625, ysize, xsize),dtype=np.complex64)
    #G = np.zeros((625, ysize, xsize),dtype=np.complex64)
    #Gprime = np.zeros((625, ysize, xsize),dtype=np.complex64)  
    
    absF = pow(diffPattern3D,0.5)
    
    # parameters
    beta = 0.9               # HIO feedback parameter
    thresholdMask1 = 0.04    # Initial threshold for mask
    thresholdMask2 = 0.2     # threshold for mask 
    sigma = 3                # initial value of gaussian width
    c = absF.shape[0]/2      # center of Gaussian what about avrundning?  #lägg till center av shape0,1,2
    
    n = 39                   #number of iterations
    
    # initialize vector for error calculation
    #err2 = np.zeros(shape=(1,n))
    #g = np.zeros((absF.shape[0], absF.shape[1], absF.shape[2]),dtype=np.complex64)
    
    # g represents the unknown object. Create initial guess for g using 
    # the autocorrelation function of the object. Start by Fourier 
    # transform the diffraction pattern:
    g = np.complex64(fft.fftshift(fft.ifftn(absF))) #finns några val som symmetric?
    # create autocorrelation function of g
    #kanske finns än funtion i python som gör det här
        # gör om till 3D. titta på normaliseringen!
    #
    g = abs(fft.fftshift(fft.ifftn(fft.fftn(g)*np.conj(fft.fftn(g)))))/(g.shape[0]*g.shape[1])
    # create a logical matrix that is the Shrinkwrap support
    support = g > thresholdMask1*g.max()
    # define iteration counter
    k = 0
    # Start of iterations 
    while k<n: 
        # every 20th iteration, update the support
        if k%20==0:   # VARFÖR SKULLE K VA NEG?
    
            # call function create2Dgaussian
            gauss3D = create3Dgaussian(sigma,c,absF.shape[0],absF.shape[1],absF.shape[2])
    
            # calculate the convolution the absolute value of the object wave 
            # with a gaussian
            # finns en djungel av olika alternativ hur man convolverar
            support = fftconvolve(gauss3D,abs(g),'same')    # ska man convolva med abs(g)?
    
            # Create a logical matrix that is the Shrinkwrap support
            support = support > thresholdMask2*support.max()                        #SÄKERT ATT DEN HITTAR 2D MAXET?      
    
            # reduce sigma with 1% every 20th iteration until a limit of 1.5
            if sigma >= 1.5:
                sigma = sigma*0.99
            
        # STEP 1: Fourier transform of g(x),  G(u) :
        G = np.complex64(fft.fftshift(fft.fftn(g)))
        
        # STEP 2: make |G| confirm with |F|
        Gprime = np.complex64(absF*np.exp(1j*np.angle(G)))
        
        # STEP 3: inverse Fourier transform      
        gprime = fft.ifftn(fft.ifftshift(Gprime))
        
        # STEP 4: See to that gprime satisfies its constraints: 
        
        # create inverse of the support (computional aid)
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
    #           for m in range(0,absF.shape[2]):
    #               if support[l,j,m] == 0:   
    #                   # update g(l,j) outside the support
    #                  
    #                   gprime[l,j,m] = g[l,j,m] - beta * gprime[l,j,m]
    #       np.disp(l)           
    #           else:
    #               # update g'(l,j) inside the support in one of 
    #               # the two following ways:
    #               # gprime[l,j,m] = g[l,j,m]
    #               gprime[l,j,m] = gprime[l,j,m]
    
        # overwrite g with result from iteration
        g = np.complex64(gprime)
        # set all negative values of g to zero
        g[g<0] = 0  
        
        #err2[0][k] = sum(sum( (abs(abs(fft.fftshift(fft.fftn(((g*np.conj(g))))))) - absF)**2)) / sum(sum(absF**2))
        np.disp(k)
        k = k+1
    
    plt.figure()    
    plt.imshow(abs(g[100,:,:]), cmap='gray')
    plt.title('reconstructed object  [100,:,:]')
    
    return g




def plotReconstruction():
    plt.figure()
    plt.subplot(221)
    plt.imshow(abs(g[:,:,100]), cmap='gray')
    plt.title('xyplane cut of 3Dcrystal z=100')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.axis('off')
    
    plt.subplot(222)
    plt.imshow(abs(g[100,:,:]), cmap='gray')
    plt.title('xz cut. y=100')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.axis('off')
    
    plt.subplot(223)
    plt.imshow(abs(g[:,100,:]), cmap='gray')
    plt.title('yz Dcrystal x=100')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.axis('off')
    
    plt.subplot(224)
    plt.imshow(abs(g[:,:,100]), cmap='gray')
    plt.title('xyplane z=zmiddle')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.axis('off')

#plotReconstruction()

def plotOriginal():
    plt.figure()
    plt.subplot(221)
    plt.imshow(crystal3D[:,:,100], cmap='gray')
    plt.title('xyplane cut of 3Dcrystal z=100')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.axis('off')
    
    plt.subplot(222)
    plt.imshow(crystal3D[100,:,:], cmap='gray')
    plt.title('xz cut. y=100')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.axis('off')
    
    plt.subplot(223)
    plt.imshow(crystal3D[:,100,:], cmap='gray')
    plt.title('yz Dcrystal x=100')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.axis('off')
    
    plt.subplot(224)
    plt.imshow(crystal3D[:,:,int(0.5*crystal3D.shape[2])], cmap='gray')
    plt.title('xyplane z=zmiddle')
    plt.xlabel(' x')
    plt.ylabel(' y')
    plt.axis('off')

#plotOriginal()


#plt.figure()
#plt.imshow(crystal3D, cmap='gray')
#plt.title('original object')



#############Test att använda bild som är skjuvad coord
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
