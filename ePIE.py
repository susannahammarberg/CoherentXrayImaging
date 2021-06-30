# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 11:48:03 2017

@author: Susanna Hammarberg
"""

import numpy as np
from numpy import fft
# Shifting property. fft of a shifted function. se s.9 avhandling Giewek...

import matplotlib.pyplot as plt
import matplotlib.animation as animation

def ePIE( n, diffSet, probe, objectFuncNy, objectFuncNx, ypixel, xpixel, positiony, positionx, Nxy ): 

    # size of probe and diffraction patterns
    ysize = diffSet.shape[1]
    xsize = diffSet.shape[2]
  
    
    # make sure it can hold complex numbers
    # intensitet ettor phase 0 är bra guess
    objectFunc = np.ones((objectFuncNy, objectFuncNx), dtype=np.complex64)
    
    objectIlluminated = np.ones(shape=(ysize, xsize),dtype=np.complex64)
    # allocating 
    g = np.zeros((ysize, xsize),dtype=np.complex64)
    gprime = np.zeros(( ysize, xsize),dtype=np.complex64)
    G = np.zeros((ysize, xsize),dtype=np.complex64)
    Gprime = np.zeros((ysize, xsize),dtype=np.complex64)   
    
    # define iteration counter for outer loop
    k = 0
    
    #figure for animation
#    fig = plt.figure()
#    plt.gca().invert_yaxis()
#    plt.ylabel(' [µm]')
#    plt.xlabel(' [µm]')
    
    # Initialize vector for animation data
#    ims = []
    
    # initialize vector for error calculation
    sse = np.zeros(shape=(n,1))
    gprimesum = 0
    # Start of ePIE iterations
    while k < n:
        # Start of inner loop: (where you iterate through all probe positions R)
        for u in range(0,Nxy):
            
            # define xposition in matrix from motorposition            
            yposition = int(np.round(positiony[u]/ypixel))    
            xposition = int(np.round(positionx[u]/xpixel))
            
            # Cut out the part of the image that is illuminated at R(=(ypos,xpos)
            objectIlluminated = objectFunc[yposition : yposition + ysize, xposition : xposition + xsize ]
                                                        
            # get the guessed wave field out from the object at position R (only at position R)
            g = objectIlluminated * probe      
        
            # fft the wave field at position R to Fourier space
            G = (fft.fftshift(fft.fft2(g)))
#            np.disp('G shape:')
#            print(G.shape)

            # make |PSI| confirm with the diffraction pattern from R
            Gprime = diffSet[u]*np.exp(1j*np.angle(G))
            
            # inverse Fourier transform  
            gprime =  ( fft.ifft2(fft.ifftshift(Gprime)))
        
            # update the total object function with the illuminated part
            objectFunc[yposition : yposition + ysize, xposition : xposition + xsize ] = objectIlluminated + (gprime-g) * np.conj(probe) / (np.max(abs(probe))**2)# probe* annars blir det att man delar med massa nollor
            
            #update probe function
#            if k%4==0:
#               probe = probe + 1 *(gprime-g) * np.conj(objectIlluminated)/ (np.max(abs(objectIlluminated))**2)
            #beta = 0.9
            #probe = probe + beta *(gprime-g) * np.conj(objectIlluminated)/ (np.max(abs(objectIlluminated))**2)
            
            ########################            
            # Further constraints:
            ########################
            
#            # constrain object amplitude to 1
#            temp_Oamp = abs(objectFunc)
#            temp_Oamp[temp_Oamp>1] = 1
#            temp = np.angle(objectFunc)
#            objectFunc = temp_Oamp * np.exp(1j* temp)
#            
#            ##constraint object phase to negative or 0
#            temp_Ophase = np.angle(objectFunc)
#            temp_Ophase[temp_Ophase>0] = 0
#            objectFunc = abs(objectFunc) * np.exp(1j* temp_Ophase)
            
            # This is for the PRTF (Absolut men du ska ju bara göra det efter att akka iteration är klara, alltså för k=n
            # Antingen gör detta eller jämför med stara G (men då skippar man ju en iteration)
#            if k==n:
#                gprimesum = gprimesum + fft.fft(fft.fft2(objectFunc*probe))

            # anim
#            im = plt.imshow(abs(objectFunc), animated=True, interpolation='none', extent=[0,6.837770297837617,0,6.825238081022181])
            ## Error estimate (sse)  Nu är alla diff mönster viktade på samma sätt. Inte så bra när de scans som är utanför provet är oviktiga/ger ingen information
            if u == int(Nxy/2):
                save_G_for_sse = abs(G)
#                sse[k] = sse[k] + sum(sum( )**2 
            # va är det här det är ju inte rätt:!:
            #sse[k] = sse[k] + sum(sum( (diffSet[k]**2 - abs(G)**2) / 65536  ))**2            

            
#        ims.append([im])
        # looking at the last scan?
        sse[k] = sum(sum( (diffSet[int(Nxy/2)]**2 - save_G_for_sse**2 ) / 65536 ))**2  #dela innanför
        k += 1        
        print('iteration ', k)                    
       
    # End of ePIE iterations
    
    # calculate PRTF:
        
    PRTF = (gprimesum/Nxy) / ( sum(diffSet) / Nxy)
    #te = np.fft.fftfreq(PRTF)
#    np.disp(PRTF)
    # define exit wave
    psi = np.zeros((Nxy,probe.shape[0],probe.shape[1]),dtype=np.complex64)    
    #transmis = np.zeros
    # iterate over all probe positions
#    for lu in range(0,nbr_scans):
#        # define xposition in matrix from motorposition            
#        yposition = int(np.round(positiony[lu]/ypixel))    
#        xposition = int(np.round(positionx[lu]/xpixel))
#                
#        # kolla om detta är rätt
        #psi[lu] = probe * objectFunc[yposition : yposition + ysize, xposition : xposition + xsize ]
         
 #  ani = animation.ArtistAnimation(fig, ims, interval=150, blit=True,repeat_delay=200)
    ani = 1
    return (objectFunc, probe, ani, sse, psi, PRTF)
