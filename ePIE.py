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

def ePIE( diffSet, probe, objectFuncNy, objectFuncNx, ypixel, xpixel, positiony, positionx, nbr_scans ): 

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
    
    # define iteration counter for outer loop
    k = 0
    # number of iterations of outer loop
    n = 3
    
    # figure for animation
    #fig = plt.figure()
    #plt.gca().invert_yaxis()
    #fig = plt.xlabel(' [µm]')
    plt.ylabel(' [µm]')
    plt.xlabel(' [µm]')
    
    # Initialize vector for animation data
   # ims = []
    
    # initialize vector for error calculation
    sse = np.zeros(shape=(n,1))
   
    # Start of ePIE iterations
    while k < n:
        # Start of inner loop: (where you iterate through all probe positions R)
        for u in range(0,nbr_scans):
            
            # define xposition in matrix from motorposition            
            yposition = int(np.round(positiony[u]/ypixel))    
            xposition = int(np.round(positionx[u]/xpixel))
            
            # Cut out the part of the image that is illuminated at R(=(ypos,xpos)
            objectIlluminated = objectFunc[yposition : yposition + ysize, xposition : xposition + xsize ]
                                                        
            # get the guessed wave field out from the object at position R (only at position R)
            g = objectIlluminated * probe      
        
            # fft the wave field at position R to Fourier space
            G = (fft.fftshift(fft.fft2(g)))
           
            # make |PSI| confirm with the diffraction pattern from R
            Gprime = diffSet[u]*np.exp(1j*np.angle(G))
            
            # inverse Fourier transform  
            gprime =  ( fft.ifft2(fft.ifftshift(Gprime)))
        
            # update the total object function with the illuminated part
            objectFunc[yposition : yposition + ysize, xposition : xposition + xsize ] = objectIlluminated + (gprime-g) * np.conj(probe) / (np.max(abs(probe))**2)# probe* annars blir det att man delar med massa nollor
            ##update probe function
#            if k%4==0:
#               probe = probe + 1 *(gprime-g) * np.conj(objectIlluminated)/ (np.max(abs(objectIlluminated))**2)
            beta = 0.5
            probe = probe + beta *(gprime-g) * np.conj(objectIlluminated)/ (np.max(abs(objectIlluminated))**2)
            
            # anim
#            im = plt.imshow(abs(objectFunc), animated=True, interpolation='none', extent=[0,6.837770297837617,0,6.825238081022181])

#            
#            plt.ylabel(' [µm]')
#            plt.title('Object phase')
            
   #         ims.append([im])

        sse[k] = sum(sum( (diffSet[nbr_scans-1]**2 - abs(G)**2 ) / 65536 ))**2  #dela innanför
        k=k+1        
        np.disp(k)                    
        
        #SSE[0][k] =  sum(sum(abs(Gprime - diffSet[3] )**2 ))

       
    # End of iterations
    # calculate PRTF
    # define exit wave
    psi = np.zeros((nbr_scans,probe.shape[0],probe.shape[1]),dtype=np.complex64)
    # iterate over all probe positions
    for u in range(0,nbr_scans):
        psi[u] = probe * objectFunc[yposition : yposition + ysize, xposition : xposition + xsize ]
 #   plt.colorbar()
  # ani = animation.ArtistAnimation(fig, ims, interval=5, blit=True,repeat_delay=2000)
    ani = 1
    return (objectFunc, probe, ani, sse, psi)