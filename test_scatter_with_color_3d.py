# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 17:44:36 2017

@author: Sanna
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np

xlin = np.linspace(0,3,4)
ylin = np.linspace(0,40,5)
zlin = np.linspace(0,80,6)

[ X, Y, Z ] = np.meshgrid(xlin,ylin,zlin)

# NU FUNKAR DET. ORDNINGEN här var fel. men det funkade ändå, mkt jobbigt. 
obj = np.zeros((5,4,6))   #funkar olika ordningar
obj[0,0,1] = 1
     # y x z         
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X,Y,Z, c= obj, marker ='o')
  
plt.xlabel(' x [um]')
plt.ylabel(' y [um]')
ax.set_zlabel('z [um]')
#plt.zlabel(' z [um]')
plt.title('Simple object in xyz FOV.')

