# -*- coding: utf-8 -*-
"""
Created on Thursday June 29 11:40 2017
On the experiment on NanoMAX

Created from a copy of run_ePIE.
Analysis of the NanoMAX 2017 data. Read in the Merlin and PilatusK100 data.
save Merlin (so far) data as sparce matrices
Masking nad choosing roi of data. 
Bright field, dark field, DPC analysis. 
TODO: center of mass
Plotting

@author: Susanna Hammarberg

"""
#from IPython import get_ipython
#get_ipython().magic('reset -sf')

#$ ipython --gui=qt

import sys   #to collect system path ( to collect function from another directory)
sys.path.insert(0, 'C:/Users/Sanna/Desktop/CXI/Shrinkwrap') #to collect 2Dgaussian
from create2Dgaussian import create2Dgaussian
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
import matplotlib.animation as animation
import itertools as it # used to iterate in 2 staes in a loop

#import sys   #to collect system path ( to collect function from another directory)
sys.path.insert(0, 'C:\Users\Sanna\Desktop\CXI\Ptychography') #to collect ePIE
from ePIE import ePIE  
#plt.close("all")

# start at scan nbr:
scan_name_int = 195 #458 # 94 fly scan    # is updated in for loop
scan_name_string = '%d' %scan_name_int   
first_scan_nbr = scan_name_int #save this nbr for plotting

# somehow write different for flyscan or stepscan

#detector_name = Merlin pil1M  pil110K

#directory = 'D:/exp20170628_Wallentin_nanomax/JWX33/NW2/stepscan_348/scan_0348_pil100k_' 

#directory = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX29A_NW1/scan_0%d_merlin_'%scan_name_int
directory = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX33_NW2/scan_0%d_merlin_'%scan_name_int

#directory_pil100K = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX29A_NW1/scan_0%d_pil100K_'%scan_name_int 
directory_pil100K = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX33_NW2/scan_0%d_pil100K_'%scan_name_int 
# OBS change the names to flyscan/stepscan
#directory = 'D:/exp20170628_Wallentin_nanomax/NW1/scan97/scan_0097_pil1m_' 

# glöm ej att ändra i loopen också
metadata_directory = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX33_NW2/JWX33_NW2.h5' 
#metadata_directory = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX29A_NW1/JWX29A_NW1.h5' 
  
# TODO:  
# is not used in the loop that reads in that data anymore
# but is still used in the plotting . should be defined adter the loop that reads in the data. but is currently used before that look to define how large the pil100K data matrix should be
nbr_rotations = 2

#(rows:) in one flyscan #S
nbr_rows = 17 # 16 för söndagsmacrot  # (including 0000)
#flyscan x-positions:
nbr_cols = 21# 101 #(including 0000) #nbr_positionsy = 31    #21#for scan49 J.W    #31 för scan 33
#nbr_positionsx = 31
#nbr_scansy = 17#21    scan 17 med 17x13?
#nbr_scansx = 13#21

# exp. parameters for conversion between detector and object plane
energy = 9.49   #keV   
wavelength = (1.23984E-9)/energy
#pixel_det = 172E-6   # Pixel ctr to ctr distance (w) [m] #Pilatus
#z_det = 4.2   # Pil100K
#For Merlin
pixel_det = 55E-6   # Pixel ctr to ctr distance (w) [m] #Merlin

# object-detector distance? JWs skript säger 0.7
z_Merlin = 1.065  # inte rätt för söndagsmacrot
#z_pil1M = 1.235
epsilon = 2.220446049250313e-16


#TODO: kontrollera att maskning är OK
def create_mask_Merlin():
    # Alex mask:
    data = np.load('C:\Users\Sanna\Desktop\NanoMAX062017\merlin_mask(1).npy')

#    # min mask
#    probe_mask = np.ones((512,512), dtype=np.int32)
#    # Find all cold? pixels (they show minus values)
#    
#    #probe_mask = sum(one_position) > -1   # finns inga på Merlin!?
#    # OBS ,  jag gjorde min mask uppochned eftersom Merlin var uppochned
#    
#    # remove too high intensity pixelson Merlin. 
#    #probe_mask[199,321] = 0 
#    probe_mask[217,301] = 0 
#    probe_mask[199,321] = 0 
#    
#    probe_mask[320,231] = 0 
#    probe_mask[461,64] = 0 
#    probe_mask[353,32] = 0 
#    probe_mask[236,44] = 0 
#    
#    probe_mask[38,280] = 0 
#    probe_mask[420,380] = 0 
#
#    # ta bort  korset på Merlin
#    # korset är på pixlarna 255 256 i båda led
#    probe_mask[255,:] = 0      
#    probe_mask[256,:] = 0   
#    probe_mask[:,255] = 0    
#    probe_mask[:,256] = 0 
    
    return data#probe_mask#data

# Choose mask: gather mask or make mask
mask_Merlin = create_mask_Merlin()
plt.figure()
plt.imshow(mask_Merlin)

# ev matrix to rad in Pil1M data
#diffSet=np.zeros((nbr_positions, 1043, 981))  # Pil1M

# Allocate memory Pil 100K
# TODO_ cant have this it is to large
diffSet_Pil100K = np.zeros((2,nbr_rows,nbr_cols, 195, 487),dtype=np.int32)  


# Allocate memory Merlin
# u dont allocate in python! (not even np?)

# load metadata to gather motorpositions in for loop
metadata = h5py.File( metadata_directory)
# create motorpostion array to store for each scan
motorpositions_gonphi=np.zeros((nbr_rotations))
# y is sampled once per line (see JW mail called 'SV: data')
#motorpositiony = np.zeros((nbr_rotations,nbr_rows))
## there is no samsx for flyscnas (there is one but i dont know what it means. Use 'adlink_buff' for x )
#motorpositionx = np.zeros((nbr_rotations))

row_Merlin = []
list_Merlin = []# [ [], [], [], [], [] ]            #make lists inside a list li. alt kan man skriva list_Merlin.append([]) i for loopen
motorpositiony = []
motorpositionx = []
#tuple_Merlin = ((nbr_rows), )    
# maby should use tuple istead of list..
# index number called rotation
rotation = 0
# vetor defining how many pixels the diffraction pattern is shifted with (for Merlin #S458-461,496-515)
vertical_shift = [-1,-1,0,0,0,  0,0,2,1,0,  1,1,1,0,-1,  -1,-1,-1,-1,0,  -1,-1,0,0,1,  1,-1,0,1,0,   2,0,0,1,1,  1,0,0,1,1,  1,2,2,2,4,  3,3,3,3,3,   3];

# read in and save data ROIs + motorpositions. mask data
for scan_number in it.chain(range(195, 196)):#488), range(496, 515)):
    
    # replace the scans that where rerun:     
    if (scan_number == 472):  
        scan_number = 518
    elif (scan_number == 487):
        scan_number = 519
     
    # define list to save all data from 1 rotation(all rows, all postitions):
    temp_list = []

    # loop over all 16 flyscans that constitute the set
    for row in range(0, nbr_rows):  
        
        # load hdf5-file
        scan = h5py.File( directory  + str('{0:04}'.format(row)) + '.hdf5','r') # read-only

        # load and store in np array
        data_Merlin =  np.array(scan.get('/entry_0000/measurement/Merlin/data' ), dtype=np.int32)

        # Flip up and down
        # flipping the diffpatterns en och en, maybe there is an easier way
        for col in range(0, nbr_cols):
            data_Merlin[col] = np.flipud(data_Merlin[col])
        

        # mask the Merlin data
 #       data_Merlin = data_Merlin * mask_Merlin
        
        # remove the kors from the middle of the merlin detector images
#        data_Merlin[:,255,:] = data_Merlin[:,254,:]
##        #data_Merlin[:,255,:] = (data_Merlin[:,254,:] + data_Merlin[:,257,:]    )/1
##        #np.disp((data_Merlin[:,254,:] + data_Merlin[:,257,:]    ))
##        #data_Merlin[:,256,:] = (data_Merlin[:,254,:]  +data_Merlin[:,257,:]      )/2
#        data_Merlin[:,256,:] = data_Merlin[:,257,:]
#        data_Merlin[:,:,255] = data_Merlin[:,:,254]    
#        data_Merlin[:,:,256] = data_Merlin[:,:,257]    
        
        # (gör det någon skillnad om jag gör maskningen och väljer roi på denna rad istället och skippar att 
        # spara det i en 'vanlig' matris?)
        # select roi on the detector   (for the last macro)
#        roi_y = vertical_shift[rotation] + 100, vertical_shift[rotation] + 300
#        data_Merlin = data_Merlin[:, roi_y[0]:roi_y[1], 100:380] #       data_Merlin = data_Merlin[:,130:250,250:360]
#        select roi (other scans)
        #data_Merlin = data_Merlin[:, 200: 480, 100:380]
        #data_Merlin = data_Merlin[:, 250: 430, 150:330]
        data_Merlin = data_Merlin[:, 70: 250, 150:330]
        # OBS OBS Need to change the roi in the mask aswell!!!!

        # save all images as sparse matrices in a list M
        one_row_sparse_Merlin = [sparse.lil_matrix(data_Merlin[i]) for i in xrange(nbr_cols)]
        
        # lägg till M i en list för varje rad
        temp_list.append(one_row_sparse_Merlin)

        def test_function_in_function():
            print('test')
        #test_function_in_function()
        
        #load and save transmission data from pil100K:
        scan = h5py.File( directory_pil100K  + str('{0:04}'.format(row)) + '.hdf5','r') # read-only
        data_pil = scan.get('/entry_0000/measurement/Pilatus/data' ) #pilatus data
        diffSet_Pil100K[rotation][row] = np.array(data_pil)

        
        #good or bad to delete inside loop? gets overwritten if not.
        #del scan, data_Merlin, data_pil
    


    # gather motor postions from metadata (h5-file). one file for each scan #S,not one file for each flyscan
    motorpositions_directory = '/entry%s' %scan_name_string  
    
    dataset_motorposition_gonphi = metadata.get(motorpositions_directory + '/measurement/gonphi')      
    dataset_motorpositiony = metadata.get(motorpositions_directory + '/measurement/samy')
    # instead of samx, you find the motorposition in flyscans from 'adlink_buff'
    dataset_motorpositionx = metadata.get(motorpositions_directory + '/measurement/AdLinkAI_buff') 
    



    # TODO: add gonphi
    motorpositions_gonphi[rotation] = np.array(dataset_motorposition_gonphi)
    
    motorpositiony.append(np.array(dataset_motorpositiony))
    
    # Load in x motorpositions from AdLinkAI_buff (contains additional zeros)
    motorpositionx_with_zeros = np.array(dataset_motorpositionx)
    # store the positions without the zeros
    motorpositionx.append(motorpositionx_with_zeros[:,0:nbr_cols])

    
    
    # save the whole diffraction-shebank (append one whole rotation )
    list_Merlin.append(temp_list)

    rotation = rotation + 1
    # update directories to gather next scan
    scan_name_string = '%d' %scan_number
    
    #directory = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX29A_NW1/scan_0%d_merlin_'%scan_number # sundaymacro
    directory = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX33_NW2/scan_0%d_merlin_'%scan_number 
        # inte så snyggt att deklarera om directory?
        
    np.disp('rotation:')
    np.disp(rotation)

del scan, data_Merlin, one_row_sparse_Merlin, temp_list, motorpositionx_with_zeros
#del data_pil
del dataset_motorpositionx, dataset_motorpositiony, dataset_motorposition_gonphi

# need to select roi for mask aswell as i am going to send it to ePIE
mask_Merlin = mask_Merlin[70: 250, 150:330]
        
fig = plt.figure()
plt.imshow(np.log10(list_Merlin[0][9][15].toarray()), cmap = 'jet', interpolation = 'none')
plt.colorbar()

## test 'unsparse' diffraction matrix
#J=[m.toarray() for m in M]
# for a single diffraction pattern:
#list_Merlin[0][0][0].toarray()
#plt.figure()
#plt.imshow(np.log10(J[9]))



def create_mask_Pil100K():
##    probe_mask = np.ones((diffSet.shape[1],diffSet.shape[2]))
    # Find all cold pixels (they show minus values)
    # my looking at the sum of all diffraction patterns from one random scan (the first)
    sumTotalDiffSet= sum(sum(diffSet_Pil100K[0]))
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
## apply PilK100 mask
diffSet_Pil100K = diffSet_Pil100K * probe_mask_Pil100K

#objectFunc, probe, ani, sse, psi, PRTF = ePIE(k, diffSet, probe, objectFuncNy, objectFuncNx, ypixel, xpixel, positiony, positionx, nbr_scans)
def bright_field_analysis(data):
    
    photons = np.zeros((nbr_rows,nbr_cols)) 
    #max_intensity = np.sum(  np.sum(data,axis=1) , axis=1).max()   # sum over rows and columns not sum over different diffPatterns
    for row in range(0,nbr_rows):
        for col in range(0,nbr_cols):
            photons[row,col] = sum(sum(data[row, col])) #/ max_intensity
                
    return photons

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
# if sats so that it plots the right detector name       
#detector = 1
#rotation_analysis = 0
#if (detector == 0):                                               
#
#    #bright_field = bright_field_analysis(list_Merlin-...........[rotation_analysis])
#    #dpc_x, dpc_y, pol_DPC_r, pol_DPC_phi = diff_phase_contrast(list_Merlin.........)
#    analyse_detector_name_string ='Merlin'    
#else: 
#    bright_field = bright_field_analysis(diffSet_Pil100K[rotation_analysis])
#    dpc_x, dpc_y, pol_DPC_r, pol_DPC_phi = diff_phase_contrast(diffSet_Pil100K[rotation_analysis])
#    analyse_detector_name_string ='Pil100K'    
#
#    
def plot_analysis():
    
    plt.figure()
    plt.imshow(bright_field, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: Bright field on %s'%((rotation_analysis+first_scan_nbr),analyse_detector_name_string))
    #plt.xlabel('Nominal motorpositions [um]')
    #plt.ylabel('Nominal motorpositions [um]')
    plt.colorbar()
    plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_transm'%(rotation_analysis+first_scan_nbr), bbox_inches='tight')

     
    plt.figure()
    plt.imshow(dpc_x, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: Horizontal DPC on %s'%((rotation_analysis+first_scan_nbr),analyse_detector_name_string))
    #plt.xlabel('Nominal motorpositions [um]')  
    plt.colorbar()
    plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_DPCx'%(rotation_analysis+first_scan_nbr), bbox_inches='tight')
    
    plt.figure()
    plt.imshow(dpc_y, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: Vertical DPC on %s'%((rotation_analysis+first_scan_nbr),analyse_detector_name_string))
    #plt.xlabel('Nominal motorpositions [um]')
    #plt.ylabel('Nominal motorpositions [um]')
    plt.colorbar()
    plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_DPCy'%(rotation_analysis+first_scan_nbr), bbox_inches='tight')
    
    plt.figure()
    plt.imshow(pol_DPC_r, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: DPC r on %s'%((rotation_analysis+first_scan_nbr),analyse_detector_name_string))
    plt.xlabel('Nominal motorpositions [um]')
    plt.colorbar()
    plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_DPCpol_r'%(rotation_analysis+first_scan_nbr), bbox_inches='tight')

    plt.figure()    
    plt.imshow(pol_DPC_phi, cmap = 'gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: DPC phi on %s'%((rotation_analysis+first_scan_nbr),analyse_detector_name_string))
    #plt.xlabel('Nominal motorpositions [um]')
    #plt.ylabel('Nominal motorpositions [um]')
    plt.colorbar()  
    plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_DPCpol_phi'%(rotation_analysis+first_scan_nbr), bbox_inches='tight')

#plot_analysis()

def COM(rot):
    # define a vector with length of the length of roi on the detector
    roix = np.linspace(1, 487,487)
    # define a vector with length of the height of roi on the detector
    roiy = np.linspace(1,195,195)
    
    X, Y = np.meshgrid(roix,roiy)
    
    for row in range(0,nbr_rows):
        for col in range(0,nbr_cols):
            sum(sum(diffSet_Pil100K[rot][row][col][roiy][roix]*X))/sum(sum(diffSet_Pil100K[rot][row][col][roiy][roix]))


def variables_scatterplot():
    #this is defined for scatter plot:
    x = np.linspace( first_scan_nbr , first_scan_nbr + nbr_rotations, nbr_rotations)
    y = np.linspace( 130, 250, 120)
    z = np.linspace( 250, 360, 110)
    X, Y, Z = np.meshgrid(y,x,z) # rörigt!
    X1, Y1 = np.meshgrid(z,y)
#variables_scatterplot()

 
# construct a factor I0 to compensate for intensity variation in the beam 
# where I0 = 1 is the largest value
def intens_norm_PilK100(rotation):
    
    I0 = np.zeros((nbr_rows, nbr_cols))
    for row in range(0,nbr_rows):
        for col in range (0, nbr_cols):
            I0[row][col] = sum(sum(diffSet_Pil100K[rotation][row][col]))
    I0 = I0/I0.max()
    return I0
I0 = intens_norm_PilK100(0)

# TODO: this can be for both Bragg and transmission. factor np.cos(2theta) blir ju ett för vinkel 0    
#def preps_phaseretival_2DBragg(rotation):
# Run ePIE for Flyscans    
#def run_ePIE_2D(rotation):
rotation = 0
# theta should be read out from gonphi for each particular scan=?
theta = 11.1*np.pi/180 #arb

# Sizes of roi of diffraction patterns (pixels used on the Merlin detector)
Ny = list_Merlin[0][0][0].shape[0]      
Nx = list_Merlin[0][0][0].shape[1]     

# size of one pixel in objectplane
xpixel = z_Merlin*wavelength/ (Nx * pixel_det)
ypixel = z_Merlin*wavelength/ (Ny * pixel_det)

# what the width of the diffraction pattern equals to in object plan (pixlar * pixelstorlekx/y)
sizeDiffObjectx = Nx * xpixel
sizeDiffObjecty = Ny * ypixel
   
# hur långt motorn rör sig i x och yled: 
# TODO: lägga in en factor *np.cos(2*theta)
motorpositiony[rotation] = np.cos(1*theta)* motorpositiony[rotation]
motorpositionx[rotation] = np.cos(1*theta)* motorpositionx[rotation]
motorWidthx = ( motorpositionx[rotation].max() - motorpositionx[rotation].min() ) * 1E-6
motorWidthy = ( motorpositiony[rotation].max() - motorpositiony[rotation].min() ) * 1E-6

# so the size of the object function should be enough to contain: (but actually it becomes a little bit larger because i have to round of to a hole pixel)
objectFuncSizeMaxy = motorWidthy + sizeDiffObjecty
objectFuncSizeMaxx = motorWidthx + sizeDiffObjectx

# so with a pixel-size of xpixel * ypixel, the obect function should be this many pixels:
# should i use ceil!? or floor?
objectFuncNy = ceil(objectFuncSizeMaxy / ypixel)
objectFuncNx = ceil(objectFuncSizeMaxx / xpixel)

# allocate memory for object function
objectFunc = np.zeros((int(objectFuncNy), int(objectFuncNx)))

####
# handle motorpositions so that they are in a format that work in the ePIE algorithm
# They should be in a vector, not a matrix, with length equal to the amount of grid positions
#####
# the y motorposition only has one datapoint per row so i copy it (in a confusing but good way)
positiony = np.ones((nbr_cols,nbr_rows))
positiony = (positiony * motorpositiony[rotation]).transpose()
positiony = positiony.ravel()
# convert the xpositions, which are in a matrix, to a long vector
positionx = motorpositionx[rotation].ravel()

# 'normalized' motorpostions converted to meters
positiony = (positiony - positiony.min() ) *1E-6
positionx = (positionx - positionx.min() ) *1E-6

#nbr of iterations
k = 150
# calculate the total number of 'grid' points
nbr_scans = nbr_cols*nbr_rows
# each scan is on a grid of 16 rows (16 flyscans) and 101 cols 
# in total 16*101=1616 diffraction patterns. Insert these in to a 
# numpy matrix like diffSet, that I used before.
index = 0
# FÖR SISTA MACROT ( så att båda 2 peakarna får plats ): diffSet = np.zeros((1616,200,280 ))
diffSet = np.zeros((nbr_scans,Ny, Nx))
for row in range(0,nbr_rows):
    for col in range(0,nbr_cols):
        # compensate for intensity variation with the factor I0
        diffSet[index] = list_Merlin[rotation][row][col].toarray() #/ I0[row][col] 
        index = index + 1
        

# probe construction
sigmay = 5# 14.1# 14.1               # initial value of gaussian height     #Scan51 2 x 2 
sigmax = 5# 10                    # initial value of gaussian width
probe = create2Dgaussian( sigmay, sigmax, diffSet.shape[1], diffSet.shape[2])

fig = plt.figure()
plt.imshow(np.log10(sum(diffSet)), cmap = 'jet', interpolation = 'none')
plt.colorbar()

# run 2d ePIE
objectFunc, probe, ani, sse, psi, PRTF = ePIE(k, diffSet, probe, objectFuncNy, objectFuncNx, ypixel, xpixel, positiony, positionx, nbr_scans, mask_Merlin)

#    return (objectFunc, probe)
# run ePIE 2D. use rotation as input
#objectFunc, probe = run_ePIE_2D(0)
plt.figure()     #, origin="lower"                         # sets the scale on axes. 
plt.imshow( np.angle(objectFunc), cmap='jet', interpolation='none', extent=[0,objectFuncNx*xpixel*1E6, 0,objectFuncNy*ypixel*1E6])
#plt.gca().invert_yaxis() 
plt.xlabel(' [$\mu m$]')
plt.ylabel(' [$\mu m$]')
plt.title('Scan %d: Object phase'%scan_name_int)
plt.colorbar()
plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_Ophase_k%d'%(scan_name_int, k), bbox_inches='tight')
   
plt.figure()                                                            # horisontalt vertikalt. xpixel * size(objectfunc[xled])
plt.imshow(abs(objectFunc), cmap='gray', interpolation='none', extent=[0,objectFuncNx*xpixel*1E6, 0, objectFuncNy*ypixel*1E6])
plt.xlabel(' [$\mu m$]')
plt.ylabel(' [$\mu m$]')
plt.title('Scan %d: Object amplitude'%scan_name_int)
plt.colorbar()
plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_Oamp_k%d'%(scan_name_int, k), bbox_inches='tight')


plt.figure()
plt.imshow((abs(probe)), cmap='gray', interpolation='none', extent=[ 0,sizeDiffObjectx*1E6, 0,sizeDiffObjecty*1E6])
plt.xlabel(' [$\mu m$]')
plt.ylabel(' [$\mu m$]')
plt.title('Scan %d: Probe amplitude'%scan_name_int)
plt.colorbar()
plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_Pamp_k%d'%(scan_name_int, k), bbox_inches='tight')

plt.figure()                                                            # horisontalt vertikalt
plt.imshow(np.angle(probe), cmap='jet', interpolation='none', extent=[ 0,sizeDiffObjectx*1E6, 0,sizeDiffObjecty*1E6])
plt.xlabel(' [$\mu m$]')
plt.ylabel(' [$\mu m$]')
plt.title('Scan %d: Probe phase'%scan_name_int)
plt.colorbar()
plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_Pphase_k%d'%(scan_name_int, k), bbox_inches='tight')  
    
#    #plt.gca().invert_yaxis() 

#plt.show()

# save the Bragg peak in file to open it in matlab
#scipy.io.savemat('C:/Users/Sanna/Desktop/NanoMAX062017/Bragg_peak_S458_.mat', mdict={'Braggpeak': one_position_roi})


# Look how the COM directly on the diffraction pattern varies with angle. not a good COM definition.
# remove?
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

def plot_summed_patterns_one_point(row,col):
    # sum the diffrtaction patterns in one point
    summed_point = 0
    for j in range(0,nbr_rotations):
        summed_point  = summed_point + list_Merlin[j][row][col].toarray()  
    plt.figure()
    plt.imshow(np.log10(summed_point), cmap = 'hot', interpolation = 'none')
    plt.axis('off')
    plt.colorbar()
    plt.title('masked diffraction patterns sum of one position in %d scans'%nbr_rotations)
#plot_summed_patterns_one_point(8,49)   #inser col and row you want to plot


#            # anim TABORT
#            im = plt.imshow(abs(objectFunc), animated=True, interpolation='none', extent=[0,6.837770297837617,0,6.825238081022181])
#            
#                    ims.append([im])
#                    ani = animation.ArtistAnimation(fig, ims, interval=1500, blit=True,repeat_delay=2000)
#                    plt.show()
# plot diffraction patterns merlin + pilK100
def plotalot():
    
    #figure for animation
    fig = plt.figure()
    # Initialize vector for animation data
    ims = []    
    
    
    for col in range(0,82,1): #film går upp till 82 Detta är tidsaxeln på filmen
        
        # plot a single image (single rotations, single row and column)
#        plt.figure()
#        plt.imshow(np.log10(list_Merlin[i][8][49].toarray()), cmap = 'hot', interpolation = 'none')
#        plt.colorbar()
#        plt.title('Scan_nbr_%d'%(first_scan_nbr+i))
#        
        # plot the sum of the Bragg peak, summed in the 3 (?) direction        
        summed_list_Merlin_1 = 0        
        summed_list_Merlin_2 = 0        
        summed_list_Merlin_3 = 0        
        #TODO FIXA dessa tre plottar: (som JWs film)
        for rot in range(0,49):
            summed_list_Merlin_1 = summed_list_Merlin_1 + (list_Merlin[rot][6][col].toarray())
        #plt.figure()    
        im = plt.imshow(np.log10(summed_list_Merlin_1), animated=True, cmap = 'jet', interpolation = 'none')
        ims.append([im])
    ani = animation.ArtistAnimation(fig, ims, interval=500, blit=True,repeat_delay=200)  
    plt.axis('off')
    plt.show()
    # save animation:
    ani.save('dynamic_images.mp4', writer="mencoder")
    # plot the first imjage and compare before and after the shifting
    plt.figure()
    plt.imshow(np.log10(list_Merlin[0][6][50].toarray()), animated=True, cmap = 'jet', interpolation = 'none')
    
# Dessa 2 ej rätt_
#    for j in range(0,101):
#        summed_list_Merlin_2 = summed_list_Merlin_2 + (list_Merlin[i][6][col].toarray())
#    plt.figure()    
#    plt.imshow(np.log10(summed_list_Merlin_2), cmap = 'jet', interpolation = 'none')
#    plt.colorbar()  
#
#    for j in range(0,101):
#        summed_list_Merlin_3 = summed_list_Merlin_3 + (list_Merlin[i][6][col].toarray())
#    plt.figure()    
#    plt.imshow(np.log10(summed_list_Merlin_3), cmap = 'jet', interpolation = 'none')
#    plt.colorbar()  

          
        # plot a single pil100K image
#        plt.figure()
#        plt.imshow(np.log10(diffSet_Pil100K[i][8][49]), cmap = 'hot', interpolation = 'none')
#        plt.colorbar()
#        plt.title('Scan_nbr_%d'%(first_scan_nbr+i))
        
#plotalot()

#pp=list_Merlin[0][0][0].toarray()     # OBS OBS jättestor!
#print getsizeof(pp)
#print getsizeof(list_Merlin[0][0][0])
#plt.figure()
#rowsum= np.sum(one_position,axis=1)    
#sumsum= np.sum(rowsum,axis=1)    
#xlinE = np.linspace(458,487,nbr_rotations)
#plt.plot(xlinE,sumsum,'+-')
#plt.title('Summed intensity as function of scan for one position')
#   

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



# TODO: redo
def stage_phaseretrival_3dBragg(rotation):    
#def calc_pixelsize_Merlin():
#def calc_pixelsize_Pil100K():
    # theta should be read out from gonphi for each particular scan=?
    theta = 15*np.pi/180 #arb
    
    # Sizes of roi of diffraction patterns (pixels used on the Merlin detector)
    Ny = list_Merlin[0][0][0].shape[0]      
    Nx = list_Merlin[0][0][0].shape[1]     
    
    # define pixel sizes in reciprocal space
    dq1 = 2*np.pi*pixel_det /(wavelength * z_Merlin)
    dq2 = 2*np.pi*pixel_det /(wavelength * z_Merlin)
    
    # one pixel on the object plane correspongs to ... m in the object plane at distance z from the detector (x and y)
    # är det rätt med Nx och Ny och inte tvärtom?
    dr1 = 2*np.pi / ( Ny * dq1 * np.cos(theta))
    dr2 = 2*np.pi / ( Nx * dq2 * np.cos(theta))
    #dr3 = 2*np.pi / ( N3 * dq3 * np.cos(theta))  

    dr1 = 2*np.pi / ( Ny * dq1 * np.cos(theta))
    dr2 = 2*np.pi / ( Nx * dq2 * np.cos(theta))
    #dr3 = 2*np.pi / ( N3 * dq3 * np.cos(theta))
    
    # create realspace pixel sizes
    # hur blir det här? dr1 beror av dr3
    dx = dr3*np.cos(theta)
    dy = dr2
    #dz = dr1 + dr3*np.sin(theta)
    
    #this is the size of the object(=FOV(ja, ptycho-FOV)), right?
    # what the width of the diffraction pattern equals to in object plan (pixlar * pixelstorlekx/y)
    sizeDiffObjectx =  Nx * dx
    sizeDiffObjecty =  Ny * dy
    
 
    # calculate how long each step is in x and y OBS kan också vara minus
    stepSizex = np.zeros((nbr_rows, nbr_cols))
    stepSizey = np.zeros((nbr_rows,1))
    
    # inte tillräcklig nogrannhet i motorpos för att det ska vara nödvändigt att göra detta?
    # bara använd den stepsize du valde till scanet? testa båda sätten!
    for i in range(0,nbr_rows):   #gör 2 loops for diffrent nbr of scans in y and x . convert from microns to meters
        
        stepSizey[i] = (motorpositiony[rotation][i+1] - motorpositiony[rotation][i]) * 1E-6
        for cols in range(0,nbr_cols):
            # TODO: lots and less time to do it
            stepSizex[i][:] = (motorpositionx[rotation][i+1] - motorpositionx[rotation][i]) * 1E-6
###############################♠♠♠♠♠♠OLDOLDOLD♠♠♠     

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
