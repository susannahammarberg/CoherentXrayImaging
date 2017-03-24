# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 11:48:03 2017

@author: Sanna
"""

import numpy as np
from numpy import fft
# Shifting property. fft of a shifted function. se s.9 avhandling Giewek...

import matplotlib.pyplot as plt
import math    #use np.log10 for np arrays
import matplotlib.animation as animation

def ePIE( diffSet, probe, objectSizey, objectSizex ): #diffSet, xsize, ysize, nbr_scans

    # size of probe and diffraction patterns
    ysize = diffSet.shape[1]
    xsize = diffSet.shape[2]
        
    # object får den nog inte heta object. 
    # make sure it can hold complex numbers
    # intensitet ettor phase 0 är bra guess
    objectFunc = np.ones((objectSizey, objectSizex), dtype=np.complex64)
    
    objectIlluminated = np.ones(shape=(ysize, xsize),dtype=np.complex64)
    # allocating 
    g = np.zeros((625, ysize, xsize),dtype=np.complex64)
    gprime = np.zeros((625, ysize, xsize),dtype=np.complex64)
    G = np.zeros((625, ysize, xsize),dtype=np.complex64)
    Gprime = np.zeros((625, ysize, xsize),dtype=np.complex64)
    
    
    probeSize = 25   #!=probe.shapeOBS måste ändra i probe def också
    stepsize = int(probeSize * 0.4)   #för 60% överlapp  # är en float. gör om till int
    
    nbryScans = math.floor((objectSizey - ysize) / stepsize)
    nbrxScans = math.floor((objectSizex - xsize) / stepsize)
    

    # define iteration counter for outer loop
    k = 0
    # number of iterations of outer loop
    n = 50
    
    # figure for animation
    fig = plt.figure()
    
    # Initialize vector for animation data
    ims = []
    
    # initialize vector for error calculation
   # sse = np.zeros(shape=(n,1))
    # idex for iterating through the diffraction patterns
    diffSetIndex = 0
    
    # Start of ePIE iterations
    while k < n:
        # Start of inner loop: (where you iterate through all probe positions R)
        # loop over xaxis and y axis
        for ypos in range(0, nbryScans+1):           #egentligen onödigt med 2 loopar. lättare att läsa
            for xpos in range(0, nbrxScans+1):
#                
#                 # Cut out the part of the image that is illuminated at R(=(ypos,xpos)
#                 objectIlluminated[0] = objectFunc[0: 3.7 um, 0 : 3.8 um]
#                                                         med omvanlidingsfaktor
#                 objectIlluminated[1] = objectFunc[0 + stepssize[0] : ans + 3.7, 
#                                                    (where stepsize = motorposition i+1 - 1 * omvandlingsfaktor)
                 objectIlluminated = objectFunc[ypos*stepsize:ypos*stepsize+ysize, xpos*stepsize:xpos*stepsize+xsize]
                                                  
                 # get the guessed wave field out from the object at position R (only at position R)
                 g = objectIlluminated * probe      
            
                 # fft the wave field at position R to Fourier space
                 G = (fft.fftshift(fft.fft2(g)))
                
                 # make |PSI| confirm with the diffraction pattern from R
                 Gprime = diffSet[diffSetIndex]*np.exp(1j*np.angle(G))
                 
                 # inverse Fourier transform  
                 gprime =  ( fft.ifft2(fft.ifftshift(Gprime)))
         
                 # update the TOTAL object function with the illuminated part
                 # det ska vara skillnaden mellan gamla gissningen (i punkten  R) och nya 
                 # conj()
                 objectFunc[ypos*stepsize:ypos*stepsize+ysize, xpos*stepsize:xpos*stepsize+xsize] = objectIlluminated + (gprime-g) * np.conj(probe) # probe* annars blir det att man delar med massa nollor
                 # update probe function
                 probe = probe + 0.2 * (gprime-g) * np.conj(objectIlluminated)/ (np.max(abs(objectIlluminated))**2)
                 
                 # anim
                 im = plt.imshow(abs(objectFunc), animated=True)
                 ims.append([im])
    
                 diffSetIndex = diffSetIndex+1
        
        k=k+1        
    
        # sse[k-1] = sum(sum( (diffSet[224]**2 - G**2 ) / 65536 ))**2  #dela innanför
        #SSE[0][k] =  sum(sum(abs(Gprime - diffSet[3] )**2 ))
        diffSetIndex = 0
       
    # End of iterations
    
    
    #absOri=abs(image)
    #absOb=abs(objectFunc)
    ani = animation.ArtistAnimation(fig, ims, interval=150, blit=True,repeat_delay=2000)
    #aniT = animation.ArtistAnimation(figT, imsT, interval=4000, blit=True,repeat_delay=2000)
    plt.show()
    
    plt.figure()
    plt.imshow(abs(objectFunc))
    #plt.figure()
    #plt.imshow(abs(image)
    
    return abs(objectFunc)
