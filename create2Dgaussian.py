# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 10:54:48 2017

@author: Susanna Hammarberg
"""
# this version of the function is updated so that the widthening functions sigma
# can be different for x and y

import numpy as np

def create2Dgaussian(sigmay, sigmax, height, width):
    gauss2D = np.zeros(shape=(height,width))
    for m in range(1,height):   #från 1 eller 0?
       for n in range(1,width):          
            gauss2D[m,n] = np.exp(-(((m-height/2)**2)/(2*sigmay**2) + ((n-width/2)**2)/(2*sigmax**2)))
            

    return gauss2D

    
#x= create2Dgaussian(1,1,1,1)