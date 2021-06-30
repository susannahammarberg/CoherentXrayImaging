# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 15:40:49 2017

@author: Susanna Hammarberg

Update to NanoMAX 2020 data
trying to reconstruct a 2d slice from a 3D bragg ptycho data set

HOW TO RUN:
* change scan_name_int to your desired scan nbr
* change directory and metadata_directory
* choose mask directory
* change nbr_scans, nbr_scansx, nbr_scansy
* change parameters: energy etc
* chose saved mask for vogt data or create mask with function create_mask that I wrote
* choose how to cut and center the diffpatterns (diffSet = diffSet[:, 31:195, 152:336])
* under 'probe construction' choose initial values for probe. do not forget the phase part
* choose number of iterations
"""

from IPython import get_ipython
get_ipython().magic('reset -sf')

import sys   #to collect system path ( to collect function from another directory)
sys.path.insert(0, 'C:/Users/Sanna/Desktop/CXI/Shrinkwrap') #to collect 2Dgaussian

# supports sparse matrices. dont know how they work with np matrises
#import scipy.sparse as sparse
import sys   #to collect system path ( to collect function from another directory)
from numpy import fft
from ePIE import ePIE  
from create2Dgaussian import create2Dgaussian
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import hdf5plugin #for 2020 data
import h5py
from math import floor
from math import ceil
from scipy import misc # to imoport image
from scipy.optimize import curve_fit as curve_fit
#from sys import getsizeof   #se hur mkt minne variabler tar upp 
import matplotlib.animation as animation
import time

from scipy import misc # to imoport image
from mpl_toolkits.mplot3d import axes3d
from scipy import sparse 

#from mayavi import mlab

import itertools as it # used to iterate in 2 staes in a loop

date_str = time.strftime("%Y%m%d") # -%H%M%S")



#%%
def create_mask_Merlin():
    # Alex mask:
    mask_dir = 'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/beamtime_folder/process/'
    file_name = 'merlin_mask_200430_8keV.h5'
    with h5py.File(mask_dir  + file_name,'r' ) as fp:
        mask = np.array(fp['mask'])
#
#for i in range(0,2):
#    a[ a== a.max()] = 0
#    
    return mask

# Choose mask: gather mask or make mask
mask_Merlin = create_mask_Merlin()


#%% 
scan = 444
directory = 'C:/Users/Sanna/NanoMAX_May2020_rawdata_selection/raw/' 


with h5py.File(directory + '000' + str(scan) +'.h5','r' ) as fp:
    data = np.array(fp['entry/measurement/merlin/frames'][:,0:300,256:-1])
    gonphi = np.array(fp['entry/snapshot/gonphi'])

    #plt.savefig(r"C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\diffraction\diffraction_position_%d\gonphi"%(frame) + \
#                ('%.2f'%(gonphi)).replace(".","_") + '_'  + 'S' + str(scan) + '_' + date_str)

#%%
data = data * mask_Merlin[0:300,256:-1]

#%%
   
 
plt.figure()
plt.imshow(np.log10(sum(data)), cmap='magma', interpolation='none')#magma
plt.title('Scan %d all diffraction. Masked'%(scan))       
plt.colorbar()

#%%


plt.close("all")

#scan_name_int = 444  
#scan_name_string = '%d' %scan_name_int   
#directory = 'D:/Nanomax/Wallentin/JWX31C_1/pilatus_scan_%s_' %scan_name_string 
#metadata_directory ='D:/Nanomax/Wallentin/JWX31C_1/JWX31C_1.h5' 
#motorpositions_directory = '/entry%s' %scan_name_string           

#scan definitions
Ny = 13
Nx = 67
Nxy = Nx * Ny

# exp. parameters for conversion between detector and object plane
energy = 10.0 #keV
wavelength = (1.23984E-9)/energy
pixel_det = 55E-6   # Pixel ctr to ctr distance (w) [m] #Merlin
z = 1
# num. parameter
epsilon = 2.220446049250313e-16


# gather motor postions
xmotor = 'npoint_buff/x'
ymotor = 'npoint_buff/y'
with h5py.File(directory + '000' + str(scan) +'.h5','r' ) as hf:
    x = hf['entry/measurement/%s'%xmotor][()]
    y = hf['entry/measurement/%s'%ymotor][()]
    x = x.flatten()
    y = y.flatten()
    


# Trim and center the diffraction patterns around max intensity
# look att the sum of all patterns:
#def trim_center(diffSet):

#TODO: Whrite a function that finds the center of the diffraction patterns.
# Perhaps not by centering round the pixel with highest intensity but instead
# centering around the circle in the patterns. (look at a line scan profile)
# btw, this is how you get the indexfrom max value of np array
#maxIntensity_vector_idx = np.argmax(sum(diffSet),axis=None)
#a, b  = np.unravel_index(maxIntensity_vector_idx, sum(diffSet).shape)
#np.disp(a)
#np.disp(b)

#def centering():
    #sum(diffSet).max 


def transmission_counter():
    index = 0# OK to use same variable name as in other places in code since it is in a function, right?
    photons = np.zeros((Ny,Nx)) 
    max_intensity = np.sum(  np.sum(data,axis=1) , axis=1).max()   # sum over rows and columns not sum over different diffPatterns
    for row in range(0,Ny):
        for col in range(0,Nx):
            photons[row,col] = sum(sum(data[index])) / max_intensity
            index = index+1
            
    return photons
#transmission = transmission_counter()

#test_circular probe function taken from https://mail.scipy.org/pipermail/numpy-discussion/2011-January/054470.html
# used for dark field filter and also possible for initial probe definition
# Not very round...
def circular_filter(ySize, xSize, outer_radius, inner_radius):
    inner_circle = np.zeros((ySize, xSize)).astype('uint8')
    outer_circle = np.zeros((ySize, xSize)).astype('uint8')
    cx, cy = int(xSize/2), int(ySize/2) # The center of circle
    # construct outer circle
    y_outer, x_outer = np.ogrid[-outer_radius: outer_radius, -outer_radius: outer_radius]
    index_outer = x_outer**2 + y_outer**2 <= outer_radius**2
    outer_circle[cy-outer_radius:cy+outer_radius, cx-outer_radius:cx+outer_radius][index_outer] = 1
    #construct inner circle
    y_inner, x_inner = np.ogrid[-inner_radius: inner_radius, -inner_radius: inner_radius]
    index_inner = x_inner**2 + y_inner**2 <= inner_radius**2
    inner_circle[cy-inner_radius:cy + inner_radius, cx- inner_radius:cx+ inner_radius][index_inner] = 1
    
    return outer_circle - inner_circle

outer_radius = 40   #35
inner_radius = 12     # 0 till 8 gör ingenting
#dark_field_filter = circular_filter(data.shape[1],data.shape[2],outer_radius,inner_radius)

def dark_field(dark_field_filter):
    index = 0# OK to use same variable name as in other places in code since it is in a function, right?
    filtered_diffSet =  dark_field_filter*data
    dark_field = np.zeros((Ny,Nx))
    #meanIntesnity = sum(sum(diffSet))/nbr_scans
    for row in range(0,Ny):
        for col in range(0,Nx):
            
            dark_field[row,col] = sum(sum(filtered_diffSet[index])) / sum(sum(data[index]))
            index = index+1
            
    return dark_field
#dark_field_image = dark_field(dark_field_filter)

def diff_phase_contrast():
    tempy = 0
    tempx = 0
    index = 0
    diff_phasey = np.zeros((Ny,Nx))
    diff_phasex = np.zeros((Ny,Nx))
    pol_DPC_r = np.zeros((Ny,Nx))
    pol_DPC_phi = np.zeros((Ny,Nx))
    rem_bkg_x = np.zeros((Ny,Nx))
    rem_bkg_y = np.zeros((Ny,Nx))
        
    for row in range(0, Ny):
        for col in range(0, Nx):
            
            for m in range(0, data.shape[1]):
                for n in range(0, data.shape[2]):
                    tempy = tempy + (m-Ny/2) * data[index, m, n] #/ (data[index, m, n]+ 2.220446049250313e-16)
                    tempx = tempx + (n-Nx/2) * data[index, m, n]
            # spara värdet på den första pixeln:
            # detta känns onödigt krävande för då måste if satsen kollas varje gång fast jag vet vilket k jag vill ha
            if index == 0:
                bkg_x = tempx
                bkg_y = tempy
            sum_diffSet = sum(sum(data[index]))
            diff_phasey[row, col] = tempy / sum_diffSet
            diff_phasex[row, col] = tempx / sum_diffSet
            rem_bkg_x[row,col] = diff_phasex[row,col] - bkg_x # 68.25
            rem_bkg_y[row,col] = diff_phasey[row,col] - bkg_y # 62.2
            # DPC in polar coordinates. r then phi:
            pol_DPC_r[row, col] = np.sqrt( (rem_bkg_x[row,col])**2 + (rem_bkg_y[row,col])**2)    
            pol_DPC_phi[row, col] = np.arctan( rem_bkg_y[row,col] / rem_bkg_x[row,col])
            tempy = 0
            tempx = 0
            index = index + 1
    #for row in range(0, nbr_scansy):
    #    for col in range(0, nbr_scansx):
                
    return diff_phasex, diff_phasey, pol_DPC_r, pol_DPC_phi

#dpc_x, dpc_y, pol_DPC_r, pol_DPC_phi = diff_phase_contrast()

def plot_analysis():
    
    plt.figure()
    plt.imshow(transmission, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: Transmission'%scan)
    plt.xlabel('Nominal motorpositions [um]')
    plt.ylabel('Nominal motorpositions [um]')
    plt.colorbar()
    ###plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_transm'%scan_name_int, bbox_inches='tight')

    plt.figure()
    plt.imshow(dark_field_image, cmap='gray')#, interpolation='none', extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: Dark field image'%scan)    #%d #%((scan_name))
    plt.xlabel('Nominal motorpositions [um]')
    plt.colorbar()
    ###plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_DF'%scan_name_int, bbox_inches='tight')
    
    
#    plt.figure()
#    plt.imshow(sum(diffSet), interpolation='none',  extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
#    plt.title('Scan %d:  Sum diffPatterns without dark-field filter '%scan_name_int)
#    plt.xlabel('Nominal motorpositions [um]')
#    plt.colorbar()
#    
#    plt.figure()
#    plt.imshow(dark_field_filter*sum(diffSet), interpolation='none', extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
#    plt.title('Scan %d: Sum diffPatterns with dark-field filter '%scan_name_int)
#    plt.xlabel('Nominal motorpositions [um]')
#    plt.ylabel('Nominal motorpositions [um]')
#    plt.colorbar()
    
            
    plt.figure()
    plt.imshow(dpc_x, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: Differential phase constrast x'%scan)
    plt.xlabel('Nominal motorpositions [um]')  
    plt.colorbar()
    ###plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_DPCx'%scan_name_int, bbox_inches='tight')
    
    plt.figure()
    plt.imshow(dpc_y, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: Differential phase constrast y'%scan)
    plt.xlabel('Nominal motorpositions [um]')
    plt.ylabel('Nominal motorpositions [um]')
    plt.colorbar()
    ###plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_DPCy'%scan_name_int, bbox_inches='tight')
    
    plt.figure()
    plt.imshow(pol_DPC_r, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: DPC r'%scan)
    plt.xlabel('Nominal motorpositions [um]')
    plt.colorbar()
    ###plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_DPCpol_r'%scan_name_int, bbox_inches='tight')

    plt.figure()    
    plt.imshow(pol_DPC_phi, cmap = 'gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: DPC phi'%scan)
    plt.xlabel('Nominal motorpositions [um]')
    plt.ylabel('Nominal motorpositions [um]')
    plt.colorbar()  
    ###plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_DPCpol_phi'%scan_name_int, bbox_inches='tight')

#plot_analysis()

def pad_diffPatterns(Nx,Ny): #Kan dessa tex heta Nx och Ny när det finns glabala parameterar som heter det?
    padded_diffPatterns = np.zeros((Nxy, Ny, Nx))
    x = (Nx - data.shape[2]) / 2
    y = (Ny - data.shape[1]) / 2 
    for i in range(0, Nxy):
        padded_diffPatterns[i, y: y + data.shape[1], x: x+ data.shape[2]] = data[i]
    
    np.disp(Nx)
    return padded_diffPatterns

#data = pad_diffPatterns(350,350)#   350

# Sizes of centred cut and padded diffraction patterns
Ny = data.shape[1]      
Nx = data.shape[2]

# factor for defining pixel sizes in object plane
yfactor = (1/Ny)*z*wavelength
xfactor = (1/Nx)*z*wavelength
 
## calculate how long each step is in x and y OBS kan också vara minus
#stepSizex = np.zeros((nbr_scansx,1))
#stepSizey = np.zeros((nbr_scansy,1))
#for i in range(0,nbr_scansx):   #gör 2 loops for diffrent nbr of scans in y and x . convert from microns to meters
#    stepSizex[i] = (motorpositionx[i+1] - motorpositionx[i]) * 1E-6
#    stepSizey[i] = (motorpositiony[i+1] - motorpositiony[i]) * 1E-6

# stepsize in [m]
dx = 50e-9 
dy = 50e-9

# probe construction
sigmay = 14# 14.1# 14.1            # initial value of gaussian height     #Scan51 2 x 2 
sigmax = 10# 10                    # initial value of gaussian width
probe = create2Dgaussian( sigmay, sigmax, data.shape[1], data.shape[2])

#phase = np.pi/4 * circular_filter(data.shape[1],data.shape[2],40,0)

phase = np.pi/4*probe
# create complex probe
probe = probe * np.exp(1j*phase)

plt.figure()
plt.imshow(abs(probe))
plt.figure()
plt.imshow(np.angle(probe))

# size of one pixel in objectplane. (blir annorlunda för att Nx och Ny är olika)
xpixel = xfactor/pixel_det
ypixel = yfactor/pixel_det

# what the width of the diffraction pattern equals to in object plan (pixlar * pixelstorlekx/y)
sizeDiffObjectx =  Nx * xpixel
sizeDiffObjecty =  Ny * ypixel

# hur långt motorn rör sig i x och yled: 
motorWidthx = ( x.max() - x.min() ) * 1E-6
motorWidthy = ( y.max() - y.min() ) * 1E-6

# so the size of the object function should be enough to contain: (but actually it becomes a little bit larger because i have to round of to a hole pixel)
objectFuncSizeMaxy = motorWidthy + sizeDiffObjecty
objectFuncSizeMaxx = motorWidthx + sizeDiffObjectx

# so with a pixel-size of xpixel * ypixel, the obect function should be this many pixels:
    # should i use ceil!? or floor?
objectFuncNy = ceil(objectFuncSizeMaxy / ypixel)
objectFuncNx = ceil(objectFuncSizeMaxx / xpixel)

# allocate memory for object function
objectFunc = np.zeros((objectFuncNy, objectFuncNx))

# 'normalized' motorpostions converted to meters
positiony = (y - y.min() ) *1E-6
positionx = (x - x.min() ) *1E-6

## ska dessa användas ? vet ej men får samma (!?) resultat
#positiony = abs(motorpositiony - motorpositiony.max() ) *1E-6
#positionx = abs(motorpositionx - motorpositionx.max() ) *1E-6
#
## mirror diffraction patterns
#diffSet = np.fliplr(diffSet)
#
#plt.figure()                                #                        x                  y
#plt.imshow(abs(probe), cmap='gray', interpolation='none', extent=[0,sizeDiffObjectx*1E6,0,sizeDiffObjecty*1E6])
#plt.xlabel(' [µm]')
#plt.ylabel(' [µm]')
#plt.title('Initial probe amplitude')
#plt.colorbar()
#
#plt.figure()                               
#plt.imshow(np.angle(probe), cmap='gray', interpolation='none', extent=[0,sizeDiffObjectx*1E6,0,sizeDiffObjecty*1E6])
#plt.xlabel(' [µm]')
#plt.ylabel(' [µm]')
#plt.title('Initial probe phase')
#plt.colorbar()

# run ePIE for k nbr of iterations
k = 1
objectFunc, probe, ani, sse, psi, PRTF = ePIE(k, data, probe, objectFuncNy, objectFuncNx, ypixel, xpixel, positiony, positionx, Nxy)
plt.show() #show animation
#plt.figure()
#plt.imshow(np.log10(abs(fft.fftshift(fft.fft2(objectFunc)))))


#TODO someting.. PRTF?
### make ePIE function return psi (exit wave) for every position of probe
# sum over all positions?
#psi = sum(objectFunc*probe)
#fft1 = fft.fft(psi)
#plt.title('')
#               
#plt.figure()
#plt.imshow(PRTF)
#
#plt.plot(fft1)#for hist, bins=100)
#plt.show()

# function creates a gaussian  with amplitude A, center of function c, and sigma
def gauss(x, A, c, sigma):
    return A*np.exp(-(x-c)**2/(2*sigma**2))

xCol = np.linspace(0,probe.shape[1]-1,probe.shape[1])*xpixel*1E6
xRow = np.linspace(0,probe.shape[0]-1,probe.shape[0])*ypixel*1E6     #fler punkter?
yFit = gauss(xCol,152,95,2)     #sumcolumns

yCol_data = abs(probe.sum(axis=0))  # horizontal
yRow_data = abs(probe.sum(axis=1))  # vertical 

# fit probe columns to gaussian #p0 initial guesses for fitting (optional ((else == 1 1 1))
poptCol, puCol = curve_fit(gauss, xCol, yCol_data, p0=[yCol_data.max(),1,1])
# fit probe rows to gaussian
poptRow, puRow = curve_fit(gauss, xRow, yRow_data , p0=[yRow_data.max(),1,1]) 

FWHM_col = 2.35482 * poptCol[2]      #(=2 * (sqrt(2*ln(2) ) )*sigma))
FWHM_row = 2.35482 * poptRow[2]

#TODO: 2d gaussian fit
#def gauss2d(xytuple, A, cx, cy, sigmax, sigmay):
#    (x,y) = xytuple    # hur funkar detta?
#    g = A*np.exp(- ((x-cx)**2 /(2*sigmax**2) + (y-cy)**2 /(2*sigmay**2)   ))
#    return g.ravel()
#
## Create x and y indices
#x = np.linspace(0,probe.shape[1]-1,probe.shape[1])*xpixel*1E6
#y = np.linspace(0,probe.shape[0]-1,probe.shape[0])*ypixel*1E6??
#x, y = np.meshgrid(x, y)
#xytuple = (x,y);
#y2sgauss = gauss2d(xytuple, abs(probe).max(), 82, 95, 1, 1 )
#
#data2d = abs(probe)
#popt2d, pu2 = curve_fit(gauss2d, xytuple, data2d.ravel(), p0=[data2d.max(), 92, 81, 1, 1] )
##
#plt.figure()
#plt.imshow(data2d)
#plt.contour(gauss2d(xytuple,  *popt2d ).reshape(probe.shape[0], probe.shape[1]))
#plt.title('2d probe with Gaussian fit')

#plt.figure()
#plt.plot(abs(probe.sum(axis=0)), 'b+:', label='data')                                    
#plt.plot(xCol, gauss(xCol, *poptCol), 'r-', label='fit')
#plt.plot(xCol, yFit, 'g', label='manual fit')
#plt.xlabel(' [µm]')
#plt.title('Probe summed over all columns')
#plt.legend()

##############################PLOTTING################

# Make a centred line in x and y direction on the diffraction patterns
#diffSet[:,:,int(diffSet.shape[2]/2)] = 1 
#diffSet[:,int(diffSet.shape[1]/2),:] = 1

###plot the trimmed and centered sum of diffPatterns
#linex = np.linspace(0,diffSet.shape[2],diffSet.shape[2]+1)
#liney = np.linspace(0,diffSet.shape[1],diffSet.shape[1]+1)
#lineyy = diffSet.shape[1]/2 * np.ones((linex.shape))
def plot():
     
    plt.figure()
    plt.imshow(np.log10(sum(data)+1), cmap='gray', interpolation='none')#, extent=[0,diffSet.shape[1]*pixel_det*1E3, 0, diffSet.shape[2]*pixel_det*1E3] )
    plt.xlabel(' y [mm]')
    plt.ylabel(' x [mm]')
    plt.title('Scan %d: log10 of all diffraction patterns summed'%scan)   
    plt.colorbar()
    ###plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_SumdiffPatt__k%d'%(scan_name_int, k), bbox_inches='tight')
    
    #def plott():  
    plt.figure()     #, origin="lower"                         # sets the scale on axes. 
    plt.imshow( np.angle(objectFunc), cmap='jet', interpolation='none')#, extent=[0,objectFuncNx*xpixel*1E6, 0,objectFuncNy*ypixel*1E6])
    #plt.gca().invert_yaxis() 
    plt.xlabel(' [µm]')
    plt.ylabel(' [µm]')
    plt.title('Scan %d: Object phase'%scan)
    plt.colorbar()
    ##plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_Ophase_k%d'%(scan_name_int, k), bbox_inches='tight')
       
    plt.figure()                                                            # horisontalt vertikalt. xpixel * size(objectfunc[xled])
    plt.imshow(np.log10(abs(objectFunc)), cmap='gray', interpolation='none')#, extent=[0,objectFuncNx*xpixel*1E6, 0, objectFuncNy*ypixel*1E6])
    plt.xlabel(' [µm]')
    plt.ylabel(' [µm]')
    plt.title('Scan %d: Object amplitude'%scan)
    plt.colorbar()
    ##plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_Oamp_k%d'%(scan, k), bbox_inches='tight')
    
    plt.figure()
    plt.imshow(abs(probe), cmap='gray', interpolation='none')#, extent=[0,sizeDiffObjectx*1E6, 0,sizeDiffObjecty*1E6])
    plt.xlabel(' [µm]')
    plt.ylabel(' [µm]')
    plt.title('Scan %d: Probe amplitude'%scan)
    plt.colorbar()
    ##plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_Pamp_k%d'%(scan, k), bbox_inches='tight')
    
    plt.figure()                                                            # horisontalt vertikalt
    plt.imshow(np.angle(probe), cmap='gray', interpolation='none')#, extent=[ 0,sizeDiffObjectx*1E6, 0,sizeDiffObjecty*1E6])
    plt.xlabel(' [µm]')
    plt.ylabel(' [µm]')
    plt.title('Scan %d: Probe phase'%scan)
    plt.colorbar()
    ##plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_Pphase_k%d'%(scan, k), bbox_inches='tight')  
    
    plt.figure()
    plt.plot(sse)
    plt.xlabel(' iterations ')
    plt.ylabel(' SSE ')
    plt.title('Scan %d: SSE looking at central position only'%scan_name_int)
    ##plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_SSE_k%d'%(scan, k), bbox_inches='tight')
    
    plot_x = np.linspace(0, data.shape[2]-1, data.shape[2])*xpixel*1E6
    plt.figure()
    plt.plot(plot_x ,abs(probe.sum(axis=0)), 'b+:', label='data')                                    
    plt.plot(plot_x, gauss(xCol, *poptCol), 'r-', label='fit')
    #plt.plot(plot_x, yFit, 'g', label='manual fit')
    plt.xlabel(' [µm]')
    plt.ylabel('Intensity')
    plt.title('Scan %d: Probe summed over all rows. FWHM: %f µm'%(scan,FWHM_col))     #horizontal line
    plt.legend()
    ##plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_probe_row_lineplot_k%d'%(scan_name_int, k), bbox_inches='tight')

    plot_y = np.linspace(0,data.shape[1]-1,data.shape[1])*ypixel*1E6
    plt.figure()
    plt.plot(plot_y, abs(probe.sum(axis=1)), 'b+:', label='data')  # vertical line
    plt.plot(plot_y, gauss(xRow, *poptRow),'r-', label='fit')
    plt.legend() 
    plt.xlabel(' [µm]')
    plt.title('Scan %d: Probe summed over all columns. FWHM: %f µm'%(scan_name_int,FWHM_row))
    ##plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_probe_col_lineplot_k%d'%(scan_name_int, k), bbox_inches='tight')

    def normalize_0_1(array):
        array = (array - array.min()) / (array.max() - array.min())
        return array
#
#    plt.figure()
#    x_line = np.linspace(motorpositionx[0], motorpositionx[-1], nbr_scansx)
#    row_line_nbr = 15
#    plt.plot(x_line, normalize_0_1(pol_DPC_r[row_line_nbr,:]),'r+-' ,label='pol_DPC_r')
#    plt.plot(x_line, normalize_0_1(dark_field_image[row_line_nbr,:]) ,'y+-', label='DF')  # detta är horisontella profilen   
#    plt.title('Scan %d: Scaled horizontal line profiles'%scan_name_int)
#    plt.xlabel('Nominal motorpositions [um]')
#    plt.legend()
#    plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_DF_lineplot_x_k%d'%(scan_name_int, k), bbox_inches='tight')
#    
#    plt.figure()
#    y_line = np.linspace(motorpositiony[0], motorpositiony[-1], nbr_scansy)
#    col_line_nbr = 15
#    plt.plot( y_line, normalize_0_1(pol_DPC_r[:, col_line_nbr]) ,'r+-', label='pol_DPC_r')
#    plt.plot( y_line, normalize_0_1(dark_field_image[:, col_line_nbr]), 'y+-', label='DF')  # detta är vertikala profilen   
#    plt.title('Scan %d: Scaled vertical line profiles'%scan_name_int)
#    plt.xlabel('Nominal motorpositions [um]')
#    plt.legend()
#    plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_DF_lineplot_y_k%d'%(scan_name_int, k), bbox_inches='tight')

    return 0
plot()

