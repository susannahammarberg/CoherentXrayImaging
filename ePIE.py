# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 11:48:03 2017

@author: Sanna
"""

import numpy as np
from numpy import fft
# Shifting property. fft of a shifted function. se s.9 avhandling Giewek...

import matplotlib.pyplot as plt
import matplotlib.animation as animation

def ePIE( diffSet, probe, objectFuncNy, objectFuncNx, ypixel, xpixel, positiony, positionx ): 

    # size of probe and diffraction patterns
    ysize = diffSet.shape[1]
    xsize = diffSet.shape[2]
  
    # object får den nog inte heta object. 
    # make sure it can hold complex numbers
    # intensitet ettor phase 0 är bra guess
    objectFunc = np.ones((objectFuncNy, objectFuncNx), dtype=np.complex64)
    
    objectIlluminated = np.ones(shape=(ysize, xsize),dtype=np.complex64)
    # allocating 
    g = np.zeros((625, ysize, xsize),dtype=np.complex64)
    gprime = np.zeros((625, ysize, xsize),dtype=np.complex64)
    G = np.zeros((625, ysize, xsize),dtype=np.complex64)
    Gprime = np.zeros((625, ysize, xsize),dtype=np.complex64)
    
    
    #flytta ut ur funktion
    nbr_scans = 960

    # define iteration counter for outer loop
    k = 0
    # number of iterations of outer loop
    n = 14
    
    # figure for animation
   # fig = plt.figure()
    
    # Initialize vector for animation data
#    ims = []
#    
    # initialize vector for error calculation
   # sse = np.zeros(shape=(n,1))

   
    # Start of ePIE iterations
    while k < n:
        # Start of inner loop: (where you iterate through all probe positions R)
        for u in range(0,nbr_scans):
            
            yposition = int(np.round(positiony[u]/ypixel))    
            xposition = int(np.round(positionx[u]/xpixel))

            # Cut out the part of the image that is illuminated at R(=(ypos,xpos)
            objectIlluminated = objectFunc[yposition : yposition + ysize, xposition : xposition + xsize ]
            
            #objectFunc[ypos*stepsize:ypos*stepsize+ysize, xpos*stepsize:xpos*stepsize+xsize]
                                             
            # get the guessed wave field out from the object at position R (only at position R)
            g = objectIlluminated * probe      
        
            # fft the wave field at position R to Fourier space
            G = (fft.fftshift(fft.fft2(g)))
           
            # make |PSI| confirm with the diffraction pattern from R
            Gprime = diffSet[u]*np.exp(1j*np.angle(G))
            
            # inverse Fourier transform  
            gprime =  ( fft.ifft2(fft.ifftshift(Gprime)))
        
            # update the TOTAL object function with the illuminated part
            # det ska vara skillnaden mellan gamla gissningen (i punkten  R) och nya 
            # conj()
            objectFunc[yposition : yposition + ysize, xposition : xposition + xsize ] = objectIlluminated + (gprime-g) * np.conj(probe) / (np.max(abs(probe))**2)# probe* annars blir det att man delar med massa nollor
            # update probe function
            probe = probe + 1 *(gprime-g) * np.conj(objectIlluminated)/ (np.max(abs(objectIlluminated))**2)
            
            # anim
#            im = plt.imshow(abs(objectFunc), animated=True)
#            ims.append([im])
            #np.disp(u)
        
        k=k+1        
        np.disp(k)                    # va!?
        # sse[k-1] = sum(sum( (diffSet[224]**2 - G**2 ) / 65536 ))**2  #dela innanför
        #SSE[0][k] =  sum(sum(abs(Gprime - diffSet[3] )**2 ))

       
    # End of iterations
    
    
    #absOri=abs(image)
    #absOb=abs(objectFunc)
#    ani = animation.ArtistAnimation(fig, ims, interval=150, blit=True,repeat_delay=2000)
#  
    plt.figure()
    plt.imshow(np.log10(abs(objectFunc)))
    plt.figure()
    plt.imshow(abs(objectFunc))
#    plt.figure()
#    plt.imshow(abs(gprime))
    return abs(objectFunc)