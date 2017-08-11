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
import matplotlib.animation as animation
import itertools as it # used to iterate in 2 staes in a loop
plt.close("all")

# start at scan nbr:
scan_name_int = 458 #458 # 94 fly scan    # is updated in for loop
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
  
# TODO:  
# is not used in the loop that reads in that data anymore
# but is still used in the plotting . should be defined adter the loop that reads in the data. but is currently used before that look to define how large the pil100K data matrix should be
nbr_rotations = 30

#(rows:)
nbr_rows = 16  # (including 0000)
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


#TODO: kontrollera att maskning är OK
def create_mask_Merlin():
    probe_mask = np.ones((512,512))
    # Find all cold? pixels (they show minus values)
    
    #probe_mask = sum(one_position) > -1   # finns inga på Merlin!?
    
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

# ev matrix to rad in Pil1M data
#diffSet=np.zeros((nbr_positions, 1043, 981))  # Pil1M

# Allocate memory Pil 100K
diffSet_Pil100K = np.zeros((nbr_rotations,nbr_rows,nbr_cols, 195, 487),dtype=np.int32)  


# Allocate memory Merlin
# u dont allocate in python! (not even np?)

# load metadata to gather motorpositions in for loop
metadata = h5py.File( metadata_directory)
# create motorpostion array to store for each scan
motorpositions_gonphi=np.zeros((nbr_rotations))
# y is sampled once per line (see JW mail called 'SV: data')
motorpositiony = np.zeros((nbr_rotations,nbr_rows))
# there is no samsx for flyscnas (there is one but i dont know what it means. Use 'adlink_buff' for x )
motorpositionx = np.zeros((nbr_rotations))

row_Merlin = []
list_Merlin = []# [ [], [], [], [], [] ]            #make lists inside a list li. alt kan man skriva list_Merlin.append([]) i for loopen
#tuple_Merlin = ((nbr_rows), )    
# maby should use tuple istead of list..

# This is not so nice. confusing:
# from 458-488 + 496-515
# TODO: add the 2 scans replacing 472 and 487
rotation = 0
#for rotation in it.chain(range(scan_name_int, 30), range(38, 57)):   # (== from 458-488 + 496-515)  
for scan_number in it.chain(range(458, 461)):#488), range(496, 515)):
    
    # replace the scans that where rerun:     
    if (scan_number == 472):  
        scan_number = 518
    elif (scan_number == 487):
        scan_number = 519
     
    # define list to save all data from 1 rotatin(all rows, all postitions):
    temp_list = []
    for row in range(0, nbr_rows):  
        
        # load hdf5-file
        scan = h5py.File( directory  + str('{0:04}'.format(row)) + '.hdf5','r') # read-only

        # load and store in np array
        data_Merlin =  np.array(scan.get('/entry_0000/measurement/Merlin/data' ), dtype=np.int32)
        
        # mask the Merlin data
        data_Merlin = data_Merlin * mask_Merlin
        
        # select roi
        data_Merlin = data_Merlin[:,100:300,100:380] #       data_Merlin = data_Merlin[:,130:250,250:360]
   
        # gör det någon skillnad om jag gör maskningen och väljer roi på denna rad istället och skippar att 
        # spara det i en 'vanlig' matris?
        # save all images as sparse matrices in a list M
        one_row_sparse_Merlin = [sparse.lil_matrix(data_Merlin[i]) for i in xrange(nbr_cols)]
        
        # lägg till M i en list för varje rad
        temp_list.append(one_row_sparse_Merlin)


#        # load and save transmission data from pil100K:
#        scan = h5py.File( directory_pil100K  + str('{0:04}'.format(row)) + '.hdf5','r') # read-only
#        data_pil = scan.get('/entry_0000/measurement/Pilatus/data' ) #pilatus data
#        diffSet_Pil100K[rotation][row] = np.array(data_pil)

        
        #good or bad to delete inside loop? gets overwritten if not.
        #del scan, data_Merlin, data_pil
    

#        ## gather motor postion
#        motorpositions_directory = '/entry%s' %scan_name_string   
#        dataset_motorposition_gonphi = metadata.get(motorpositions_directory + '/measurement/gonphi')
#        motorpositions_gonphi[rotation] = np.array(dataset_motorposition_gonphi)
#        
#        dataset_motorpositiony = metadata.get(motorpositions_directory + '/measurement/samy')
#        # instead of samx, you find the motorposition in flyscans frmo 'adlink_buff'
#        dataset_motorpositionx = metadata.get(motorpositions_directory + '/measurement/adlink_buff') 
#        
#        motorpositiony[rotation,:] = np.array(dataset_motorpositiony) 
#        motorpositionx[rotation] = np.array(dataset_motorpositionx) 
    
    # save the whole shebank (append one whole rotation )
    list_Merlin.append(temp_list)

    rotation = rotation + 1
    # update directories to gather next scan
    scan_name_string = '%d' %scan_number
    directory = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX29A_NW1/scan_0%d_merlin_'%scan_number 
        # inte så snyggt att deklarera om directory?
        
    np.disp('rotation:')
    np.disp(rotation)

   
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
# apply PilK100 mask
diffSet_Pil100K = diffSet_Pil100K * probe_mask_Pil100K


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
#detector = 0
#rotation_analysis = 1
#if (detector == 0):                                               
#
#    #bright_field = bright_field_analysis(list_Merlin-...........[rotation_analysis])
#    #dpc_x, dpc_y, pol_DPC_r, pol_DPC_phi = diff_phase_contrast(list_Merlin.........)
#    analyse_detector_name_string ='Merlin'    
#else: 
#    bright_field = bright_field_analysis(diffSet_Pil100K[rotation_analysis])
#    dpc_x, dpc_y, pol_DPC_r, pol_DPC_phi = diff_phase_contrast(diffSet_Pil100K[rotation_analysis])
#    analyse_detector_name_string ='Pil100K'    

    
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



#this is defined for scatter plot:
x = np.linspace( first_scan_nbr , first_scan_nbr + nbr_rotations, nbr_rotations)
y = np.linspace( 130, 250, 120)
z = np.linspace( 250, 360, 110)
X, Y, Z = np.meshgrid(y,x,z) # rörigt!
X1, Y1 = np.meshgrid(z,y)


# save the Bragg peak in file to open it in matlab
#scipy.io.savemat('C:/Users/Sanna/Desktop/NanoMAX062017/Bragg_peak_S458_.mat', mdict={'Braggpeak': one_position_roi})

def shift_scans():
    
    ##scan_no = [458,459,460,461,462, 463,464,465,466,467,  468,469,470,471,518,  473,474,475,476,477,   478,479,480,481,482,  483,484,485,486,519,  488, 496,497,498, 499,   500,501,502,503,504,  505,506,507,508,509,  510,511,512,513,514,  515];  
    #% Using this reduces the FOV!
    # Vector for vertical shift of an entire scan (in pixels)
    vertical_shift = [-1,-1,0,0,0,  0,0,2,1,0,  1,1,1,0,-1,  -1,-1,-1,-1,0,  -1,-1,0,0,1,  1,-1,0,1,0,   2,0,0,1,1,  1,0,0,1,1,  1,2,2,2,4,  3,3,3,3,3,   3];
    
    for rot in range(0,1):
        list_Merlin[rot][row][col] 






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
    
    
    for col in range(0,82,1): #film går upp till 82  #(0,nbr_rotations,1) Detta är tidsaxeln på filmen
        
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
        
plotalot()


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


#def calc_pixelsize_Merlin():
## Sizes of roi of diffraction patterns
#Ny = list_Merlin[0][0][0].shape[0]      
#Nx = list_Merlin[0][0][0].shape[1]     
#    
##def calc_pixelsize_Pil100K():
##
### Sizes of roi of diffraction patterns
##Ny = list_Merlin[0][0][0].shape[0]      
##Nx = list_Merlin[0][0][0].shape[1]      
#
## factor for defining pixel sizes in object plane
#yfactor = (1./Ny)*z*wavelength
#xfactor = (1./Nx)*z*wavelength
# 
## calculate how long each step is in x and y OBS kan också vara minus
#stepSizex = np.zeros((nbr_scansx,1))
#stepSizey = np.zeros((nbr_scansy,1))
#for i in range(0,nbr_scansx):   #gör 2 loops for diffrent nbr of scans in y and x . convert from microns to meters
#    stepSizex[i] = (motorpositionx[i+1] - motorpositionx[i]) * 1E-6
#    stepSizey[i] = (motorpositiony[i+1] - motorpositiony[i]) * 1E-6


## size of one pixel in objectplane. (blir annorlunda för att Nx och Ny är olika)
#xpixel = xfactor/pixel_det
#ypixel = yfactor/pixel_det
#
## what the width of the diffraction pattern equals to in object plan (pixlar * pixelstorlekx/y)
#sizeDiffObjectx =  Nx * xpixel
#sizeDiffObjecty =  Ny * ypixel
#
## hur långt motorn rör sig i x och yled: 
#motorWidthx = ( motorpositionx.max() - motorpositionx.min() ) * 1E-6
#motorWidthy = ( motorpositiony.max() - motorpositiony.min() ) * 1E-6
#
## so the size of the object function should be enough to contain: (but actually it becomes a little bit larger because i have to round of to a hole pixel)
#objectFuncSizeMaxy = motorWidthy + sizeDiffObjecty
#objectFuncSizeMaxx = motorWidthx + sizeDiffObjectx
#
## so with a pixel-size of xpixel * ypixel, the obect function should be this many pixels:
#    # should i use ceil!? or floor?
#objectFuncNy = ceil(objectFuncSizeMaxy / ypixel)
#objectFuncNx = ceil(objectFuncSizeMaxx / xpixel)
#
## allocate memory for object function
#objectFunc = np.zeros((int(objectFuncNy), int(objectFuncNx)))
#
## 'normalized' motorpostions converted to meters
#positiony = (motorpositiony - motorpositiony.min() ) *1E-6
#positionx = (motorpositionx - motorpositionx.min() ) *1E-6