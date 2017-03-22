# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 11:48:03 2017
@author: Sanna
"""
import numpy as np
from numpy import fft
# Shifting property. fft of a shifted function. se s.9 avhandling Giewek...
from scipy import misc
import matplotlib.pyplot as plt
import math    #use np.log10 for np arrays
import matplotlib.animation as animation


#image = misc.imread('gubbe.gif')
#image = image[:,:,0]
#image = misc.imread('greens.jpg',flatten=True) #flatten gives a single gray scale layer
image = misc.imread('fruit.jpg',flatten=True) #flatten gives a single gray scale layer
#image = misc.imread('ko400400.jpg',flatten=True) #flatten gives a single gray scale layer
#image = image[:,:,0]

# generalisea ':
image = np.zeros(shape=(400, 400),dtype=np.complex64) +  (0.2*image/image.max()) +  (0.8*1j*image/image.max())  # bara complex del nu

Nsize = 256  #(brukade vara 60)

imageSet=np.zeros((625,Nsize, Nsize),dtype=np.complex64) #have to uese double paranthesis
diffSet=np.zeros((625, Nsize, Nsize),dtype=np.complex64) #have to uese double paranthesis
# diffSet bilder kommer ju inte var komplexa från experiment

#minnesallokering
# object får den nog inte heta object. 
# make sure it can hold complex numbers
# intensitet ettor phase 0 är bra guess
objectFunc = np.ones(image.shape,dtype=np.complex64)
objectIlluminated = np.ones(shape=(Nsize, Nsize),dtype=np.complex64)

g = np.zeros((625,Nsize, Nsize),dtype=np.complex64)
gprime = np.zeros((625,Nsize, Nsize),dtype=np.complex64)
G = np.zeros((625,Nsize, Nsize),dtype=np.complex64)
Gprime = np.zeros((625,Nsize, Nsize),dtype=np.complex64)

probeSize = 25   #!=probe.shapeOBS måste ändra i probe def också
stepsize = int(probeSize * 0.4)   #för 60% överlapp  # är en float. gör om till int
nbrScans = math.floor((image.shape[0] - Nsize) / stepsize)

# define probe for diffraction patterns
probe = np.zeros(shape=(Nsize, Nsize),dtype=np.complex64)
probe[116:141, 116:141] = 1  #np.ones(shape=(20, 20)). OK


# klipp ut submatriserna till diffraktionsmönstreerna
# and create diffraction patterns (sqrt of diff patterns) by fft2 of O*P
index=0
for ypos in range(0, nbrScans+1):   #+1?
    for xpos in range(0, nbrScans+1): # +1         
        
        imageSet[index] = image[ypos*stepsize: ypos*stepsize+Nsize, xpos*stepsize: xpos*stepsize+Nsize]
        diffSet[index] =  abs(fft.fftshift(fft.fft2(imageSet[index]*probe)))
        index = index+1
# look att diff patterns:        
#plt.imshow(np.log10(abs(diffSet[5])))

#imageSet = None # or tot delete it 
del imageSet


# update probe to initial guess. size and value
probe[116:141, 116:141] = 0.1  #np.ones(shape=(20, 20)). OK
probe[110:151, 110:151] = 0.1  #np.ones(shape=(20, 20)). OK

# define iteration counter for outer loop
k = 0
# number of iterations of outer loop
n = 15

#fig = plt.figure()
#figT = plt.figure()
ims = []
#imsT = []

# initialize vector for error calculation
sse = np.zeros(shape=(n,1))
error = np.zeros(shape=(n,1))
diffSetIndex = 0
while k < n:
    # Start of inner loop: (where you iterate through all probe positions R)
    # loop over xaxis and y axis
    for ypos in range(0, nbrScans+1):           #egentligen onödigt med 2 loopar. lättare att läsa
        for xpos in range(0, nbrScans+1):
            
             # Cut out the part of the image that is illuminated at R(=(ypos,xpos)
             objectIlluminated = objectFunc[ypos*stepsize:ypos*stepsize+Nsize, xpos*stepsize:xpos*stepsize+Nsize]
                                                  
             # get the guessed wave field out from the object at position R (only at position R)
             g = objectIlluminated * probe      
        
             # fft the wave field at position R to Fourier space
             G = (fft.fftshift(fft.fft2(g)))
            
             # make |PSI| confirm with the diffraction pattern from R
             Gprime = diffSet[diffSetIndex]*np.exp(1j*np.angle(G))
             
             # inverse Fourier transform  
                       #fft.ifftshift? igen?
             gprime = fft.ifft2(fft.ifftshift(Gprime))
     
             # update the TOTAL object function with the illuminated part
             # det ska vara skillnaden mellan gamla gissningen (i punkten  R) och nya 
             # if a function other than just ones are used for the probe it should be
             # adressed here (divide by maxvalue)
             objectFunc[ypos*stepsize:ypos*stepsize+Nsize, xpos*stepsize:xpos*stepsize+Nsize] =  objectIlluminated + (gprime-g) * np.conj(probe) # probe* annars blir det att man delar med massa nollor
             
             # update probe function
             probe =  probe + 0.2*(gprime-g) * np.conj(objectIlluminated)/ (np.max(abs(objectIlluminated))**2) 
             
             # anim
#             im = plt.imshow(abs(objectFunc), animated=True)
#             ims.append([im])

             diffSetIndex = diffSetIndex+1

    error[k] = sum(sum((gprime-g).real))
  #  sse[k] = sum(sum( (diffSet[224]**2 - G**2 ))) / 65536 2  #dela innanför
    k=k+1   
     
    # sse[k-1] = sum(sum( (diffSet[224]**2 - G**2 ) / 65536 ))**2  #dela innanför
    #SSE[0][k] =  sum(sum(abs(Gprime - diffSet[3] )**2 ))
    diffSetIndex = 0
# End of iterations
   
#absOri=abs(image)
#absOb=abs(objectFunc)

#ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True,repeat_delay=2000)
#plt.show()

plt.figure()
plt.imshow(abs(objectFunc))
plt.figure()
plt.imshow(np.log10(abs(probe)))
#plt.figure()
#plt.imshow(abs(image)