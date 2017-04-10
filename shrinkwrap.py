# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:04:38 2017

@author: HonkyT
"""
#from numpy import fft # om man skriver såhär begöver man itne använda np.
# men eftersom jag behöver fler saker från nupy och jag redan skrivit np på alla gör jag inte det.
import numpy as np
from scipy import signal
from create2Dgaussian import create2Dgaussian  #import all functions

data = np.genfromtxt("diffractionPattern.txt", dtype='float')
#fattar inte hur man skriver det som function på ett smidigt sätt.
# går ju inte att kolla på variablerna om man kör som funktion
#def Shrinkwrap(diffPattern):
    # the intensity of the diffraction pattern is F^2
    
absF = pow(data,0.5)

# parameters
beta = 0.9               # HIO feedback parameter
thresholdMask1 = 0.04;   # Initial threshold for mask
thresholdMask2 = 0.2;    # threshold for mask 
sigma = 3;               # initial value of gaussian width
c = absF.shape[0]/2 #varför blir det en float?      # center of Gaussian

n=420 #number of iterations

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
    
#    # STEP 4: See to that gprime satisfies its constraints: 
    for l in range(1,absF.shape[0]):
       for j in range(1,absF.shape[1]):         ## ska DET VARA RANGE 0 ÄR FÖRSTA INDEXEDT 0?
           if support[l,j] == 0:   
               # update g(l,j) outside the support
               gprime[l,j] = g[l,j] - beta * gprime[l,j]
               
           else:
               # update g'(l,j) inside the support in one of 
               # the two following ways:
               #gprime(l,j) = g(l,j);
               gprime[l,j] = gprime[l,j]

    # overwrite g with result from iteration
    g = gprime
    # set all negative values of g to zero
    g[g<0] = 0  

    # för att lätt kunna titta på g
    realg = g.real
    
    err2[0][k] = sum(sum( (abs(abs(np.fft.fftshift(np.fft.fft2(((g*np.conj(g))))))) - absF)**2)) / sum(sum(absF**2))

    k = k+1