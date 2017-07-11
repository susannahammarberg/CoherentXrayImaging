# -*- coding: utf-8 -*-
"""
Created on Thursday June 29 11:40 2017
On the experiment on NanoMAX

Created from a copy of run_ePIE.
Trying to read in data from the Bragg rocking curve
and plot a 3D peak
@author: Susanna Hammarberg

"""
#from IPython import get_ipython
#get_ipython().magic('reset -sf')

#$ ipython --gui=qt

#import sys   #to collect system path ( to collect function from another directory)
#sys.path.insert(0, 'C:/Users/HonkyT/Desktop/CXI/Shrinkwrap') #to collect 2Dgaussian

# supports sparse matrices. dont know how they work with np matrises
#import scipy.sparse as sparse
#%matplotlib qt5
from numpy import fft
import matplotlib.pyplot as plt
import numpy as np
import h5py
from math import floor
from math import ceil
import scipy
from scipy import misc # to imoport image
from mpl_toolkits.mplot3d import axes3d
#from mayavi import mlab

plt.close("all")

# start at scan nbr:
scan_name_int = 458 # 94 fly scan    # is updated in for loop
scan_name_string = '%d' %scan_name_int   
first_scan_nbr = scan_name_int #save this nbr for plotting

# somehow write different for flyscan or stepscan

#detector_name = Merlin pil1M  pil110K

#directory = 'D:/exp20170628_Wallentin_nanomax/JWX33/NW2/stepscan_348/scan_0348_pil100k_' 

directory = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX29A_NW1/scan_0%d_merlin_'%scan_name_int
directory_pil100K = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX29A_NW1/scan_0%d_pil100K_'%scan_name_int 
# OBS change the names to flyscan/stepscan
#directory = 'D:/exp20170628_Wallentin_nanomax/NW1/scan97/scan_0097_pil1m_' 

#metadata_directory = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX33/NW2/JWX33_NW2.h5' 
metadata_directory = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX29A_NW1/JWX29A_NW1.h5' 

      

nbr_rotations = 5

#(rows:)
nbr_rows = 16  # (including 0000)
#flyscan x-positions:
nbr_cols = 101 #(including 0000) #nbr_positionsy = 31    #21#for scan49 J.W    #31 för scan 33
#nbr_positionsx = 31
#nbr_scansy = 17#21    scan 17 med 17x13?
#nbr_scansx = 13#21

# exp. parameters for conversion between detector and object plane
energy = 9.49   
wavelength = (1.23984E-9)/energy
pixel_det = 172E-6   # Pixel ctr to ctr distance (w) [m] #Pilatus
z_det = 4.2   # Pil100K
#z_Merlin = 1.065  
#z_pil1M = 1.235
epsilon = 2.220446049250313e-16


# create matrix to hold raw diffraction patterns
#diffSet=np.zeros((nbr_positions, 1043, 981))  # Pil1M

diffSet_one_row = np.zeros((nbr_cols, 195, 487))   # Pil100K
diffSet_one_position= np.zeros((nbr_rotations, 195, 487))#,dtype=np.complex64)

one_row = np.zeros(( nbr_cols, 512, 512))#,dtype=np.complex64)   # Merlin
one_position= np.zeros((nbr_rotations, 512, 512))#,dtype=np.complex64)

# load metadata to gather motorpositions in for loop
metadata = h5py.File( metadata_directory)
# create motorpostion array to store for each scan
motorpositions_gonphi=np.zeros((nbr_rotations))
motorpositiony = np.zeros((nbr_rotations,nbr_rows))
motorpositionx = np.zeros((nbr_rotations))

# read data from hdf5-files
for rotation in range(0, nbr_rotations):    # could equallt be scan_nbr instead of 'rotation'
    
#    # to load all rows in a flyscan:
#    for position in range(0,nbr_positions): 
    
#        # what is the name of pil1M if pil100K is called Pilatus?
#        data_position = scan.get('/entry_0000/measurement/Merlin/data' )
#        #diffSet = np.array(data_scan)  
#        diffSet[position] = np.array(data_position)          
#    del scan, data_position    

    # to load a single row in a flyscan:
    row = 8
    position = row
    
    scan = h5py.File( directory  + str('{0:04}'.format(position)) + '.hdf5','r') # read-only
    #data_scan = scan.get('/entry_0000/measurement/Pilatus/data' ) #pilatus data
    # what is the name of pil1M if pil100K is called Pilatus?
    data_scan = scan.get('/entry_0000/measurement/Merlin/data' )
    one_row = np.array(data_scan)
    
    # select one position from the row
    col = 49
    one_position[rotation] = one_row[col,:,:]       
    
    # load and save transmission data from pil100K:
    scan = h5py.File( directory_pil100K  + str('{0:04}'.format(position)) + '.hdf5','r') # read-only
    data_position = scan.get('/entry_0000/measurement/Pilatus/data' ) #pilatus data
    diffSet_one_row = np.array(data_position)
    
    diffSet_one_position[rotation] = diffSet_one_row[col,:,:]
    
    del scan, data_scan
    ## gather motor postion
    motorpositions_directory = '/entry%s' %scan_name_string   
    dataset_motorposition_gonphi = metadata.get(motorpositions_directory + '/measurement/gonphi')
    motorpositions_gonphi[rotation] = np.array(dataset_motorposition_gonphi)
    
    dataset_motorpositiony = metadata.get(motorpositions_directory + '/measurement/samy')
    dataset_motorpositionx = metadata.get(motorpositions_directory + '/measurement/samx')
    
    motorpositiony[rotation,:] = np.array(dataset_motorpositiony) 
    motorpositionx[rotation] = np.array(dataset_motorpositionx) 
    
    # update directories to gather next scan
    scan_name_int = scan_name_int + 1
    scan_name_string = '%d' %scan_name_int
    directory = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX29A_NW1/scan_0%d_merlin_'%scan_name_int 
    # inte så snyggt att deklarera om detta?
    
    np.disp('rotation:')
    np.disp(rotation)
# delete inside loop? why? no
#del one_row

def create_mask_Merlin():
##    probe_mask = np.ones((diffSet.shape[1],diffSet.shape[2]))
    # Find all cold? pixels (they show minus values)
    
    probe_mask = sum(one_position) > -1   # finns inga på Merlin!?
    
    # remove too high intensity pixelson Merlin. 
    #probe_mask[199,321] = 0 
    probe_mask[217,301] = 0 
    probe_mask[199,321] = 0 
    
    probe_mask[320,231] = 0 
    probe_mask[461,64] = 0 
    probe_mask[353,32] = 0 
    probe_mask[236,44] = 0 
    
    probe_mask[38,280] = 0 
    probe_mask[420,380] = 0 
        
    return probe_mask

# Choose mask: gather mask or make mask
mask_Merlin = create_mask_Merlin()
# apply mask
one_position = one_position * mask_Merlin


# def ROI
one_position_roi = one_position[:,130:250,250:360]

#this is definedfor the scatter plot:
x = np.linspace( first_scan_nbr , first_scan_nbr + nbr_rotations, nbr_rotations)
y = np.linspace( 130, 250, 120)
z = np.linspace( 250, 360, 110)
X, Y, Z = np.meshgrid(y,x,z) # rörigt!
X1, Y1 = np.meshgrid(z,y)


# save the Bragg peak in file to open it in matlab
scipy.io.savemat('C:/Users/Sanna/Desktop/NanoMAX062017/Bragg_peak_S458_.mat', mdict={'Braggpeak': one_position_roi})


#plt.figure()
#plt.imshow(sum(one_position_roi))


plt.figure()
plt.imshow(np.log10((np.sum(one_position,axis=0))), cmap = 'hot', interpolation = 'none')
plt.axis('off')
plt.colorbar()
plt.title('masked diffraction patterns sum of one position in 30 scans')

for i in range(0,nbr_rotations,1):
    plt.figure()
    plt.imshow(np.log10(one_position[i]), cmap = 'hot', interpolation = 'none')
    plt.colorbar()
    plt.title('Scan_nbr_%d'%(first_scan_nbr+i))
    plt.figure()
    plt.imshow(np.log10(diffSet_onepoint[i]), cmap = 'hot', interpolation = 'none')
    plt.colorbar()
    plt.title('Scan_nbr_%d'%(first_scan_nbr+i))
    
    
#plt.figure()
#rowsum= np.sum(one_position,axis=1)    
#sumsum= np.sum(rowsum,axis=1)    
#xlinE = np.linspace(458,487,30)
#plt.plot(xlinE,sumsum,'+-')
#plt.title('Summed intensity as function of scan for one position')
#   
#plt.figure()
#plt.imshow((one_position[20]), cmap = 'hot', interpolation = 'none')
#plt.title('Scan 478')
#plt.colorbar()
##plt.title('Scan_nbr_%d'%(first_scan_nbr+i))
#
#plt.figure()
#plt.imshow(np.log10(one_position[21]), cmap = 'hot', interpolation = 'none')
#plt.title('Scan 479')
#plt.colorbar()
#
#plt.figure()
#plt.imshow((one_position[22]), cmap = 'hot', interpolation = 'none')
#plt.title('Scan 480')
#plt.colorbar()


# TODO: gör en 3D scatter plot med färkodning för intensiteten:
    # gör meshgrids? lnispaces för xyz andvänd one_position för färg
def plot_Bragg():
#    #scatter 3 ritar ut en ring för varje punkt specificerad av vektorerna (x,y,z)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #cmap = one
    ax.scatter(X, Y, Z, cmap ='jet' , c = np.log10(one_position_roi), marker ='o', alpha =0.05)
    #rätt!???:
    plt.ylabel(' rotation')
    plt.xlabel(' yDet')
    ax.set_zlabel('zDet ')
    #plt.zlabel(' zDet')
    plt.title('Bragg peak')
#plot_Bragg() 


def funkar_ej():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #X, Y, Z = axes3d.get_test_data(0.05)
    #cset = ax.contour(X, Y, one_position[1])
    #ax.contourf(X1,Y1, .1*np.sin(3*X)*np.sin(5*Y))
    
    levels = np.linspace(-1, 1, 5)
    
    ax.contourf(Y1,X1, one_position_roi[0], cmap ='hot', levels=.1*levels, zdir='z')
    ax.contourf(Y1,X1,20+ one_position_roi[8], cmap ='hot', levels=3+.1*levels, zdir='z')
    ax.contourf(Y1,X1, one_position_roi[16], cmap ='hot', levels=5+.1*levels, zdir='z')
    #ax.clabel(cset, fontsize=9, inline=1)
    
    plt.show()



#plt.title('Scan_nbr_%d'%(first_scan_nbr+i))
#plt.figure()
#plt.imshow(np.log10((np.sum(diffSet[:,0,:,:], axis=0))), cmap = 'hot', interpolation = 'none')
#plt.axis('off')
#plt.colorbar()

#plt.imshow(np.log10(sum(sum(diffSet)))) # for the flyscans you sum sum
##plt.figure()
##scan_nbr= 110
##plt.plot(np.log10(np.sum(diffSet[scan_nbr],axis=0)))
##plt.title('Scan_nbr_%d'%scan_nbr)
##
### Eg there are 201 steps on the rocking curve
###rocking_steps= np.linspace(1,200)
##summed_diffSet = np.zeros((201))
##for i in range(0,201):
##
##    summed_diffSet[i] = sum(sum(diffSet[i]))
##
##    
##plt.figure()
##plt.plot(motorpositions_gonphi, np.log10(summed_diffSet))
##plt.figure()
##plt.imshow(np.log10(sum(diffSet)), interpolation='none')
###D:\exp20170628_Wallentin_nanomax
##arr = diffSet
##scipy.io.savemat('C:/Users/Sanna/Desktop/NanoMAX062017/tt.mat', mdict={'arr': arr})
##
### maximum intensity at at 168.08