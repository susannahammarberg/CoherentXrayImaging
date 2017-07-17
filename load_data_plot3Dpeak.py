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
import matplotlib.colors as colors
import numpy as np
import h5py
from math import floor
from math import ceil
#import scipy
from scipy import misc # to imoport image
from mpl_toolkits.mplot3d import axes3d
from scipy import sparse 
from sys import getsizeof   #se hur mkt minne variabler tar upp  
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

     

nbr_rotations = 2

#(rows:)
nbr_rows = 8  # (including 0000)
#flyscan x-positions:
nbr_cols = 101 #(including 0000) #nbr_positionsy = 31    #21#for scan49 J.W    #31 för scan 33
#nbr_positionsx = 31
#nbr_scansy = 17#21    scan 17 med 17x13?
#nbr_scansx = 13#21

# exp. parameters for conversion between detector and object plane
energy = 9.49   #keV   
wavelength = (1.23984E-9)/energy
pixel_det = 172E-6   # Pixel ctr to ctr distance (w) [m] #Pilatus
z_det = 4.2   # Pil100K
#z_Merlin = 1.065  
#z_pil1M = 1.235
epsilon = 2.220446049250313e-16


# create matrix to hold raw diffraction patterns
# one hole Merlin set:

#getsizeof(nbr_rows) # returns values in bytes
#getsizeof(diffSet_Merlin)  

#diffSet=np.zeros((nbr_positions, 1043, 981))  # Pil1M

# Allocate memory Pil 100K
diffSet_Pil100K = np.zeros((nbr_rows,nbr_cols, 195, 487),dtype=np.int32)  
diffSet_one_row = np.zeros((nbr_cols, 195, 487))   # Pil100K
diffSet_one_position= np.zeros((nbr_rotations, 195, 487))#,dtype=np.complex64)
# Allocate memory Merlin

diffSet_Merlin = np.zeros((nbr_rows,nbr_cols, 512, 512),dtype=np.int32)  
one_row = np.zeros(( nbr_cols, 512, 512))#,dtype=np.complex64)   # Merlin
one_position= np.zeros((nbr_rotations, 512, 512))#,dtype=np.complex64)

# load metadata to gather motorpositions in for loop
metadata = h5py.File( metadata_directory)
# create motorpostion array to store for each scan
motorpositions_gonphi=np.zeros((nbr_rotations))
motorpositiony = np.zeros((nbr_rotations,nbr_rows))
motorpositionx = np.zeros((nbr_rotations))

# insert all diffraction patterns for all positions,  rows and rotations into one tuple
#diffSet_Merlin = []
#diffSet_Merlin = sparse.csc_matrix((30, nbr_rows,nbr_cols, 512, 512))
list_Merlin= [ [], [], [], [], [] ]            #make lists inside a list li
#tuple_Merlin = ((nbr_rows), )    
# maby should use tuple istead of list..

for rotation in range(0, nbr_rotations):    # could equallt be scan_nbr instead of 'rotation'
    # load all rows 0 : nbr_rows:
    for row in range(0, nbr_rows): 
        position = row        
    #        # what is the name of pil1M if pil100K is called Pilatus?
    
        scan = h5py.File( directory  + str('{0:04}'.format(position)) + '.hdf5','r') # read-only
        #data_scan = scan.get('/entry_0000/measurement/Pilatus/data' ) #pilatus data
        # what is the name of pil1M if pil100K is called Pilatus?
        data_scan = scan.get('/entry_0000/measurement/Merlin/data' )
        # save a hole set (all postions) diffraction patterns but for only one angle
        diffSet_Merlin[row] = np.array(data_scan)
        
        one_row = np.array(data_scan)
        
        # select one position from the row
        col = 35 #49
              
        
        one_position[rotation] = one_row[col,:,:]       
        sparse_Merlin = sparse.csc_matrix(one_position[0])
        # OBS: random: insert one rotations on position from one row into a list
        
        # ett värde fölr varje rad (ej alla positioner är med : )
        #TODO:
        #tuple_Merlin[row] = ((sparse_Merlin) )
        # load and save transmission data from pil100K:
        scan = h5py.File( directory_pil100K  + str('{0:04}'.format(position)) + '.hdf5','r') # read-only
        data_pil = scan.get('/entry_0000/measurement/Pilatus/data' ) #pilatus data
        diffSet_Pil100K[row] = np.array(data_pil)
        diffSet_one_row = np.array(data_pil)
        
        diffSet_one_position[rotation] = diffSet_one_row[col,:,:]
        
        del scan, data_scan, data_pil
#        ## gather motor postion
#        motorpositions_directory = '/entry%s' %scan_name_string   
#        dataset_motorposition_gonphi = metadata.get(motorpositions_directory + '/measurement/gonphi')
#        motorpositions_gonphi[rotation] = np.array(dataset_motorposition_gonphi)
#        
#        dataset_motorpositiony = metadata.get(motorpositions_directory + '/measurement/samy')
#        dataset_motorpositionx = metadata.get(motorpositions_directory + '/measurement/samx')
#        # instead of samx, you find the motorposition in flyscans frmo 'adlink_buff'
#        
#        motorpositiony[rotation,:] = np.array(dataset_motorpositiony) 
#        motorpositionx[rotation] = np.array(dataset_motorpositionx) 
        
        # update directories to gather next scan

    scan_name_int = scan_name_int + 1
    scan_name_string = '%d' %scan_name_int
    directory = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX29A_NW1/scan_0%d_merlin_'%scan_name_int 
        # inte så snyggt att deklarera om detta?
    list_Merlin[0].append(sparse_Merlin)    
    np.disp('rotation:')
    np.disp(rotation)
# delete inside loop? why? no
#del one_row
#plt.figure()
#plt.imshow(np.log10(np.sum(np.sum(diffSet_Merlin[:,:,:,:],axis=0),axis=0)), interpolation='none')
#plt.colorbar()
#plt.title('Scan_nbr_%d'%(scan_name_int))
#
#plt.figure()
#plt.subplot(221)
#plt.imshow(np.log10(np.sum(diffSet_Merlin[:,3,:,:],axis=0)), interpolation='none')
#
#plt.subplot(222)
#plt.imshow(abs(fft.fftshift(fft.fft2(((np.sum(diffSet_Merlin[:,3,:,:],axis=0)))))), interpolation='none')

#TODO: kontrollera att maskning är OK
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
diffSet_Merlin = diffSet_Merlin * mask_Merlin

def create_mask_Pil100K():
##    probe_mask = np.ones((diffSet.shape[1],diffSet.shape[2]))
    # Find all cold pixels (they show minus values)
    sumTotalDiffSet= sum(sum(diffSet_Pil100K))
    probe_mask = sumTotalDiffSet > 0
    
    # remove too high intensity pixels. 
    probe_mask[84:116, 225:258] = 0 
    probe_mask[152,229] = 0 
#    j=239
#    probe_mask[111,j:245] = 0 
#    probe_mask[112,j:245] = 0 
#    probe_mask[113,j:245] = 0  
#    probe_mask[114,j:245] = 0 
#    probe_mask[115,j:245] = 0 
#    probe_mask[113,242] = 0 #pixel with highest intensity
    return probe_mask

probe_mask_Pil100K = create_mask_Pil100K()
# apply PilK100 mask
diffSet_Pil100K = diffSet_Pil100K * probe_mask_Pil100K
#plt.figure()
#plt.imshow((sum(diffSet_Pil100K[10])))


def bright_field(data):
    
    photons = np.zeros((nbr_rows,nbr_cols)) 
    #max_intensity = np.sum(  np.sum(diffSet_Merlin,axis=1) , axis=1).max()   # sum over rows and columns not sum over different diffPatterns
    for row in range(0,nbr_rows):
        for col in range(0,nbr_cols):
            photons[row,col] = sum(sum(data[row, col])) #/ max_intensity
                
    return photons
#choose what detector to do bright field on
bright_field = bright_field(diffSet_Pil100K)

def diff_phase_contrast(data):
    tempy = 0
    tempx = 0
    
    diff_phasey = np.zeros((nbr_rows,nbr_cols))
    diff_phasex = np.zeros((nbr_rows,nbr_cols))
    pol_DPC_r = np.zeros((nbr_rows,nbr_cols))
    pol_DPC_phi = np.zeros((nbr_rows,nbr_cols))
    rem_bkg_x = np.zeros((nbr_rows,nbr_cols))
    rem_bkg_y = np.zeros((nbr_rows,nbr_cols))
        
    for row in range(0, nbr_rows):
        for col in range(0, nbr_cols):
            
            # m and n pixels in the diffraction data
            for m in range(0, data.shape[2]):
                for n in range(0, data.shape[3]):
                    tempy = tempy + (m-nbr_rows/2) * data[row, col, m, n] #/ (diffSet[index, m, n]+ 2.220446049250313e-16)
                    tempx = tempx + (n-nbr_cols/2) * data[row, col, m, n]
            # spara värdet på den första pixeln:
            # detta känns onödigt krävande för då måste if satsen kollas varje gång fast jag vet vilket k jag vill ha
            if row == 0 and col == 0:
                bkg_x = tempx
                bkg_y = tempy
            sum_diffSet = sum(sum(data[row,col]))
            diff_phasey[row, col] = tempy / sum_diffSet
            diff_phasex[row, col] = tempx / sum_diffSet
            rem_bkg_x[row,col] = diff_phasex[row,col] - bkg_x # 68.25
            rem_bkg_y[row,col] = diff_phasey[row,col] - bkg_y # 62.2
            # DPC in polar coordinates. r then phi:
            pol_DPC_r[row, col] = np.sqrt( (rem_bkg_x[row,col])**2 + (rem_bkg_y[row,col])**2)    
            pol_DPC_phi[row, col] = np.arctan( rem_bkg_y[row,col] / rem_bkg_x[row,col])
            tempy = 0
            tempx = 0
        print row     
    return diff_phasex, diff_phasey, pol_DPC_r, pol_DPC_phi

# Choose which detector images to analyze:
        # TODO: fix this with the string name
analyse_detector_name_string ='Merlin'
#dpc_x, dpc_y, pol_DPC_r, pol_DPC_phi = diff_phase_contrast(diffSet_%s)%analyse_detector_name_string
#dpc_x, dpc_y, pol_DPC_r, pol_DPC_phi = diff_phase_contrast(diffSet_Pil100K)

def plot_analysis():
    
    plt.figure()
    plt.imshow(bright_field, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: Bright field'%scan_name_int)
    #plt.xlabel('Nominal motorpositions [um]')
    #plt.ylabel('Nominal motorpositions [um]')
    plt.colorbar()
    #plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_transm'%scan_name_int, bbox_inches='tight')

     
    plt.figure()
    plt.imshow(dpc_x, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: Horizontal DPC'%scan_name_int)
    #plt.xlabel('Nominal motorpositions [um]')  
    plt.colorbar()
    #plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_DPCx'%scan_name_int, bbox_inches='tight')
    
    plt.figure()
    plt.imshow(dpc_y, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: Vertical DPC'%scan_name_int)
    #plt.xlabel('Nominal motorpositions [um]')
    #plt.ylabel('Nominal motorpositions [um]')
    plt.colorbar()
    #plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_DPCy'%scan_name_int, bbox_inches='tight')
    
    plt.figure()
    plt.imshow(pol_DPC_r, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: DPC r'%scan_name_int)
    plt.xlabel('Nominal motorpositions [um]')
    plt.colorbar()
    #plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_DPCpol_r'%scan_name_int, bbox_inches='tight')

    plt.figure()    
    plt.imshow(pol_DPC_phi, cmap = 'gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: DPC phi'%scan_name_int)
    #plt.xlabel('Nominal motorpositions [um]')
    #plt.ylabel('Nominal motorpositions [um]')
    plt.colorbar()  
    #plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_DPCpol_phi'%scan_name_int, bbox_inches='tight')

#plot_analysis()

# def ROI
one_position_roi = one_position[:,130:250,250:360]

# select an even smaller roi:
#test=one_position_roi[15,54:64,47:57]


#plt.figure()
#plt.subplot(221)
#plt.imshow(one_position_roi[15])
#plt.subplot(222)
#plt.imshow(abs(fft.fftshift(fft.ifft2(fft.ifftshift(one_position_roi[15])))))

#this is defined to be able to plot scatter plot:
x = np.linspace( first_scan_nbr , first_scan_nbr + nbr_rotations, nbr_rotations)
y = np.linspace( 130, 250, 120)
z = np.linspace( 250, 360, 110)
X, Y, Z = np.meshgrid(y,x,z) # rörigt!
X1, Y1 = np.meshgrid(z,y)


# save the Bragg peak in file to open it in matlab
#scipy.io.savemat('C:/Users/Sanna/Desktop/NanoMAX062017/Bragg_peak_S458_.mat', mdict={'Braggpeak': one_position_roi})


#plt.figure()
#plt.imshow(sum(one_position_roi))


plt.figure()
plt.imshow(np.log10((np.sum(one_position,axis=0))), cmap = 'hot', interpolation = 'none')
plt.axis('off')
plt.colorbar()
plt.title('masked diffraction patterns sum of one position in %d scans'%nbr_rotations)

#Look how the COM directly on the diffraction pattern varies with angle. not a good COM definition.
def COM_variation(j, nbr_iter):

    for i in range (j,nbr_iter):
        xindex = np.argmax(np.sum(one_position[i],axis=0))
        yindex = np.argmax(np.sum(one_position[i],axis=1))
        reddot=np.zeros((512,512))
            
        # Make a centred line in x and y intersection at COM
        reddot[:,xindex] = 500000 
        reddot[yindex,:] = 500000 
        np.disp( xindex)
        plt.figure()
        noes  = ['spring', 'autumn']
        plt.imshow(np.log10(one_position[i]), cmap=noes[1] , interpolation = 'none')
        plt.imshow(np.log10(reddot))
        #plt.imshow(np.log10(one_position[1]), cmap = 'hot', interpolation = 'none')
        #plt.colorbar() funkar ej med flera imshows
        plt.title('Scan_nbr_%d'%(first_scan_nbr+i))
        
#COM_variation(0,3)    

# plot diffraction patterns merlin + pilK100
#for i in range(0,nbr_rotations,1):   #(0,nbr_rotations,1)
#    plt.figure()
#    plt.imshow(np.log10(one_position[i]), cmap = 'hot', interpolation = 'none')
#    plt.colorbar()
#    plt.title('Scan_nbr_%d'%(first_scan_nbr+i))
#    plt.figure()
#    plt.imshow(np.log10(diffSet_one_position[i]), cmap = 'hot', interpolation = 'none')
#    plt.colorbar()
#    plt.title('Scan_nbr_%d'%(first_scan_nbr+i))
    
    
plt.figure()
rowsum= np.sum(one_position,axis=1)    
sumsum= np.sum(rowsum,axis=1)    
xlinE = np.linspace(458,487,nbr_rotations)
plt.plot(xlinE,sumsum,'+-')
plt.title('Summed intensity as function of scan for one position')
   

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
#summed_diffSet = np.zeros((201))
#for i in range(0,201):
#
#    summed_diffSet[i] = sum(sum(diffSet[i]))
#
#    
##plt.figure()
##plt.plot(motorpositions_gonphi, np.log10(summed_diffSet))
##plt.figure()
##plt.imshow(np.log10(sum(diffSet)), interpolation='none')
###D:\exp20170628_Wallentin_nanomax
##arr = diffSet
##scipy.io.savemat('C:/Users/Sanna/Desktop/NanoMAX062017/tt.mat', mdict={'arr': arr})
##
### maximum intensity at at 168.08