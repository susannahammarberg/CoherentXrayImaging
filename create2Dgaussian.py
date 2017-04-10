# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 10:54:48 2017

@author: Susanna Hammarberg
"""
import numpy as np

def create2Dgaussian(sigma, center, width, length):
    gauss2D = np.zeros(shape=(width,length))
    for m in range(0,width): 
       for n in range(0,length):          
            gauss2D[m,n] = np.exp(-(((n-center)**2)+(m-center)**2)/(2*sigma**2))
            

    return gauss2D

    
#x= create2Dgaussian(1,1,1,1)