# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 10:54:48 2017

@author: Susanna Hammarberg
"""
import numpy as np

def create3Dgaussian(sigma, center, width, length, height):
    gauss3D = np.zeros(shape=(width,length,height))
    for m in range(0,width): 
       for n in range(0,length):
           for o in range(0,height):
               # missing normalization. good?
               gauss3D[m,n,o] = np.exp(-(((n-center)**2)+(m-center)**2+(o-center)**2)/(2*sigma**3)) 
            

    return gauss3D

    
#x= create2Dgaussian(1,1,1,1)