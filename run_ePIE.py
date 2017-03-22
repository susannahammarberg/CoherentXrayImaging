# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 10:24:25 2017

@author: HonkyT
"""
from scipy import misc
from ePIE import ePIE   
import numpy as np
from numpy import fft

from math import floor    #use np.log10 for np arrays  used

#image = misc.imread('gubbe.gif', flatten=True)
#image = misc.imread('greens.jpg',flatten=True)
image = misc.imread('fruit.jpg',flatten=True)
#image = misc.imread('ko400400.jpg',flatten=True) #flatten gives a single gray scale layer
#image = image[:,:,0]
objectSizey = image.shape[0]
objectSizex = image.shape[1]

# make a complex and a real part of the image
image = np.zeros((objectSizey, objectSizex), dtype=np.complex64) + (0.2*image/image.max()) + (0.8*1j*image/image.max())

nbrScans = 625#441# 625
# size of probe and diffraction patterns
probeSizey = 256#50#256
probeSizex = 256#50# 256

probeInnerSize = 25   #!=probe.shapeOBS måste ändra i probe def också
stepsize = int(probeInnerSize * 0.4)   #för 60% överlapp  # är en float. gör om till int

nbryScans = floor((objectSizey - probeSizey) / stepsize)
nbrxScans = floor((objectSizex - probeSizex) / stepsize)
 
# probe must be the same size as the object for multiplication
# initial guess for the probe size and value
probe = np.zeros(shape=(probeSizey, probeSizex),dtype=np.complex64)
#probe[116:141, 116:141] = 1  #np.ones(shape=(20, 20)). OK
probe[floor(probeSizey/2 - probeInnerSize/2) :floor(probeSizey/2 - probeInnerSize/2)+probeInnerSize,  floor(probeSizex/2 - probeInnerSize/2):floor(probeSizex/2 - probeInnerSize/2)+probeInnerSize] = 1


# create the diffraction patterns 
imageSet=np.zeros((nbrScans, probeSizey, probeSizex),dtype=np.complex64) #have to uese double paranthesis
diffSet=np.zeros((nbrScans, probeSizey , probeSizex),dtype=np.complex64) #have to uese double paranthesis
 
# diffSet bilder kommer ju inte var komplexa från experiment
    
index=0
for ypos in range(0, nbryScans+1):   #+1? Japp
    for xpos in range(0, nbrxScans+1): # +1         
        
        imageSet[index] = image[ypos*stepsize: ypos*stepsize+ probeSizey, xpos*stepsize: xpos*stepsize+probeSizex]
        diffSet[index] =  abs(fft.fftshift(fft.fft2(imageSet[index]*probe)))
        index = index+1
# look att diff patterns:        
#plt.imshow(np.log10(abs(diffSet[5])))
#imageSet = None # or tot delete it 
del imageSet
    
    
    
#animation = ePIE(diffSet, probe, objectSizey, objectSizex, probeSizey, probeSizex)




