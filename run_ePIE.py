# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 15:40:49 2017

@author: Susanna Hammarberg

See s.50 Giewekemeyer thesis för att gå mellan detektorplan och objektplan
"""
#from IPython import get_ipython
#get_ipython().magic('reset -sf')

import sys   #to collect system path ( to collect function from another directory)
sys.path.insert(0, 'C:/Users/HonkyT/Desktop/CXI/Shrinkwrap') #to collect 2Dgaussian

# supports sparse matrices. dont know how they work with np matrises
#import scipy.sparse as sparse
from numpy import fft
from ePIE import ePIE  
from create2Dgaussian import create2Dgaussian
import matplotlib.pyplot as plt
import numpy as np
import h5py
from math import floor
from math import ceil
from scipy import misc # to imoport image
from scipy.optimize import curve_fit as curve_fit
#from sys import getsizeof   #se hur mkt minne variabler tar upp 
#import matplotlib.animation as animation

plt.close("all")

scan_name = 51
#directory = 'C:/Users/HonkyT/Desktop/CXI/Ptychography/scan33/pilatus_scan_33_'
#directory = 'F:/Nanomax/Wallentin/JWX31C_1/pilatus_scan_38_'   #stående nanotråd
#directory = 'F:/Nanomax/Wallentin/JWX31C_1/pilatus_scan_17_'
directory = 'F:/Nanomax/Wallentin/JWX31C_1/pilatus_scan_51_'
#directory = 'F:/Nanomax/Wallentin/JWX31C_1/pilatus_scan_49_'
#directory = 'F:/Nanomax/Vogt_ptycho/scan33/pilatus_scan_33_'
#directory = 'F:/Nanomax/Vogt_ptycho/scan83/pilatus_scan_83_'

#metadata_directory = 'F:/Nanomax/Vogt_ptycho/scan83/DiWCr4_2.h5' #2 metadatafiler, den första är 1kB
#metadata_directory = 'F:/Nanomax/Vogt_ptycho/scan33/DiWCr4_1.h5' 
metadata_directory ='F:/Nanomax/Wallentin/JWX31C_1/JWX31C_1.h5' 

# No mask created for Wallentin (by alex)
#mask_directory = 'scan33_mask.hdf5'
#mask_directory = 'scan83_mask.hdf5'

motorpositions_directory = '/entry51' #J.W  
#motorpositions_directory = '/entry38' #J.W
#motorpositions_directory = '/entry17' #J.W
#motorpositions_directory = '/entry49' #J.W  
#motorpositions_directory = '/entry83'  
#motorpositions_directory = '/entry33'  


nbr_scans = 961#961#441  #221 scan17    #441     # 441 Scan49,38 J.W    #961 scan 33
nbr_scansy = 31    #21#for scan49 J.W    #31 för scan 33
nbr_scansx = 31
#nbr_scansy = 17#21    scan 17 med 17x13?
#nbr_scansx = 13#21

# create matrix to hold raw diffraction patterns
diffSet=np.zeros((nbr_scans, 195, 487))   

# read data from hdf5-files
for scan_nbr in range(0,nbr_scans): 
    scan = h5py.File( directory  + str('{0:04}'.format(scan_nbr)) + '.hdf5','r') # read-only
    data_scan = scan.get('/entry_0000/measurement/Pilatus/data' )
    diffSet[scan_nbr] = np.array(data_scan)   #Varför har den tre dimensioner?

del scan, data_scan

# TODO look at differences between diff patterns. remove an average value from each image.

def create_mask():
##    probe_mask = np.ones((diffSet.shape[1],diffSet.shape[2]))
    # Find all cold? pixels (they show minus values)
    sumTotalDiffSet= sum(diffSet)
    probe_mask = sumTotalDiffSet > 0
    
    #probe_mask[6,59] = 0
    #probe_mask[6,60] = 0
    #probe_mask[2,243] = 0
    #probe_mask[2,244] = 0
    #probe_mask[11,304] = 0
    #probe_mask[11,305] = 0
    #probe_mask[56,101] = 0
    
    # remove too high intensity pixels
#    probe_mask[111,241] = 0
#    probe_mask[112,240] = 0
#this
    j=239

    probe_mask[111,j:245] = 0
    probe_mask[112,j:245] = 0
    probe_mask[113,j:245] = 0 
    probe_mask[114,j:245] = 0
    probe_mask[115,j:245] = 0
    #tothis
#    
#    probe_mask[112,241] = 0
#    probe_mask[112,242] = 0
#    probe_mask[112,243] = 0
#    probe_mask[113,241] = 0 
    probe_mask[113,242] = 0 #pixel with highest intensity
#    probe_mask[113,243] = 0
#    probe_mask[114,241] = 0
#    probe_mask[114,242] = 0
#    probe_mask[114,243] = 0

    return probe_mask

# Choose mask: gather mask or make mask'
#
#mask_file = h5py.File(mask_directory)
#meta_mask = mask_file.get('/mask')
#probe_mask = np.array(meta_mask)

probe_mask = create_mask()

# apply mask
diffSet = probe_mask * diffSet

del probe_mask

plt.figure()
plt.imshow(np.log10(sum(diffSet)), cmap='gray', interpolation='none')
plt.title('log10 diffSet sum uncut')
plt.colorbar()

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

# inte riktigt rätt för Wallentin, max är vid y=82 (rätt) x=92 (hä, 96)
#diffSet = diffSet[:, 31:195, 150:342] # Vogt33
diffSet = diffSet[:, 31:195, 152:336] # Wallentin51                ##242
#diffSet = diffSet[:, 31:195, 152:330] # Wallentin38
#diffSet = diffSet[:, 31:195, 146:336]  # Wallentin17 rätt?

Ny = diffSet.shape[1]        # nbr_pixels of centred and cut diffraction patterns
Nx = diffSet.shape[2]

def transmission_counter():
    index = 0# OK to use same variable name as in other places in code since it is in a function, right?
    photons = np.zeros((nbr_scansy,nbr_scansx)) 
    max_intensity = np.sum(  np.sum(diffSet,axis=1) , axis=1).max()   # sum over rows and columns not sum over different diffPatterns
    for row in range(0,nbr_scansy):
        for col in range(0,nbr_scansx):
            photons[row,col] = sum(sum(diffSet[index])) / max_intensity
            index = index+1
            
    return photons
transmission = transmission_counter()

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

outer_radius = 35
inner_radius = 16
dark_field_filter = circular_filter(diffSet.shape[1],diffSet.shape[2],outer_radius,inner_radius)

def dark_field(dark_field_filter):
    index = 0# OK to use same variable name as in other places in code since it is in a function, right?
    filtered_diffSet =  dark_field_filter*diffSet
    dark_field = np.zeros((nbr_scansy,nbr_scansx))
    #meanIntesnity = sum(sum(diffSet))/nbr_scans
    for row in range(0,nbr_scansy):
        for col in range(0,nbr_scansx):
            
            dark_field[row,col] = sum(sum(filtered_diffSet[index])) / sum(sum(diffSet[index]))
            index = index+1
            
    return dark_field

def diff_phase_contrast():
    tempy = 0
    tempx = 0
    index = 0
    diff_phasey = np.zeros((nbr_scansy, nbr_scansx))
    diff_phasex = np.zeros((nbr_scansy, nbr_scansx))
    pol_DPC = np.zeros((nbr_scansy,nbr_scansx))
    test_x = np.zeros((nbr_scansy,nbr_scansx))
    test_y = np.zeros((nbr_scansy,nbr_scansx))
        
    for row in range(0, nbr_scansy):
        for col in range(0, nbr_scansx):
            
            for m in range(0, diffSet.shape[1]):
                for n in range(0, diffSet.shape[2]):
                    tempy = tempy + (m-nbr_scansy/2) * diffSet[index, m, n] #/ (diffSet[index, m, n]+ 2.220446049250313e-16)
                    tempx = tempx + (n-nbr_scansx/2) * diffSet[index, m, n]
            # spara värdet på den första pixeln:
            # detta känns onödigt krävande för då måste if satsen kollas varje gång fast jag vet vilket k jag vill ha
#            if index == 0:
#                bkg_x = tempx
#                bkg_y = tempy
            # ska jag dela med sum sum?
            diff_phasey[row, col] = tempy / sum(sum(diffSet[index]))
            diff_phasex[row, col] = tempx / sum(sum(diffSet[index]))
            test_x[row,col] = (tempx ) / sum(sum(diffSet[index])) - 68.25
            test_y[row,col] = (tempy ) / sum(sum(diffSet[index])) - 62.2
            # DPC in polar coordinates. r:
            pol_DPC[row, col] = np.sqrt(((tempy ) / sum(sum(diffSet[index])) - 62.2 )**2 + ((tempx ) / sum(sum(diffSet[index])) - 68.25)**2)
            # pol_DPC[row, col] = np.sqrt((tempy)**2 ... gjorde jag innan men det är ju fel då glömmer jag ju att dela. minns ej varför jag gör det men det är nog det som är rätt
            tempy = 0
            tempx = 0
            index = index + 1
            
    plt.figure()    
    plt.imshow(test_x, cmap = 'gray', interpolation='none')
    plt.title('DPCx bkg compensated')
    plt.colorbar()   
    plt.figure()    
    plt.imshow(test_y, cmap = 'gray', interpolation='none')
    plt.title('DPCy bkg compensated')
    plt.colorbar()          
    return diff_phasex, diff_phasey, pol_DPC

dpc_x, dpc_y, pol_DPC = diff_phase_contrast()

dark_field_image = dark_field(dark_field_filter)

# TODO: check if correct
#def cart2pol(x,y):    
#    r = np.sqrt(x**2 + y**2)
#    phi = np.arctan2(y,x)
#    return(r, phi)

#rCord, phiCord = cart2pol() #[1,2,3], [2,4,6]
    

def plot_analysis():
    
    plt.figure()
    plt.imshow(transmission, cmap='gray', interpolation='none')
    plt.title('Transmission')
    plt.colorbar()

    plt.figure()
    plt.imshow(dark_field_image, cmap='gray', interpolation='none')
    plt.title('Scan %d: Dark field image'%((scan_name)))
    plt.colorbar()
    
    plt.figure()
    plt.imshow(sum(diffSet), interpolation='none')
    plt.title('Sum of all diffraction patterns (centred and cut) without a dark-field filter ')
    plt.colorbar()
    
    plt.figure()
    plt.imshow(dark_field_filter*sum(diffSet), interpolation='none')
    plt.title('Sum of all diffraction patterns (centred and cut) with a dark-field filter ')
    plt.colorbar()
            
    plt.figure()
    plt.imshow(dpc_x, cmap='gray', interpolation='none')
    plt.title('Differential phase constrast x')
    plt.colorbar()
    
    plt.figure()
    plt.imshow(dpc_y, cmap='gray', interpolation='none')
    plt.title('Differential phase constrast y')
    plt.colorbar()
    plt.show()
    
    plt.figure()
    plt.imshow(pol_DPC, cmap='gray', interpolation='none')
    plt.title('Differential phase constrast in pol cord r')
    plt.colorbar()
    plt.show()

plot_analysis()

# TODO: write function for padding of diff Patterns med storled som inparameterar Nx Ny
def pad_diffPatterns(Nx,Ny): #Kan dessa tex heta Nx och Ny när det finns glabala parameterar som heter det?
    np.disp(Nx)
    return 0

pad_diffPatterns(4,5)
    
    

# parameters for conversion between detector and object plane
energy = 10.72#     # ?Wallentin.    
#energy = 8.17 #vogt # keV  
wavelength = (1.23984E-9)/energy
pixel_det = 172E-6   # Pixel ctr to ctr distance (w) [m] #Pilatus
z = 4.3

# factor for defining pixel sizes in object plane
yfactor = (1/Ny)*z*wavelength
xfactor = (1/Nx)*z*wavelength

# gather motor postions
metadata = h5py.File( metadata_directory)
dataset_motorpositiony = metadata.get(motorpositions_directory + '/measurement/samy')
dataset_motorpositionx = metadata.get(motorpositions_directory + '/measurement/samx')
motorpositiony = np.array(dataset_motorpositiony) 
motorpositionx = np.array(dataset_motorpositionx) 

del metadata
 
# calculate how long each step is in x and y OBS kan också vara minus
stepSizex = np.zeros((nbr_scansx,1))
stepSizey = np.zeros((nbr_scansy,1))
for i in range(0,nbr_scansx):   #gör 2 loops for diffrent nbr of scans in y and x . convert from microns to meters
    stepSizex[i] = (motorpositionx[i+1] - motorpositionx[i]) * 1E-6
    stepSizey[i] = (motorpositiony[i+1] - motorpositiony[i]) * 1E-6

# probe construction
sigmay = 1 #2 14.1# 14.1               # initial value of gaussian height
sigmax = 1 #2                    # initial value of gaussian width
probe = create2Dgaussian( sigmay, sigmax, diffSet.shape[1], diffSet.shape[2])

phase = np.pi/4 * circular_filter(diffSet.shape[1],diffSet.shape[2],1,0)

#phase = np.zeros((diffSet.shape[1],diffSet.shape[2]))
 ##initial guess for probe (square with amplitude 1 everywhere):
#probeInnerSize = 12#with ones
#probe = np.zeros((diffSet.shape[1],diffSet.shape[2]), dtype=np.complex64)
#probe[116:141, 116:141] = 1  #np.ones(shape=(20, 20)). OK
## kolla om detta blev rätt
#probe[floor(diffSet.shape[1]/2 - probeInnerSize/2) :floor(diffSet.shape[1]/2 - probeInnerSize/2)+probeInnerSize,  floor(diffSet.shape[2]/2 - probeInnerSize/2):floor(diffSet.shape[2]/2 - probeInnerSize/2)+probeInnerSize] = 1

#probeInnerSize = 12#with ones
#phase = np.zeros((diffSet.shape[1],diffSet.shape[2]), dtype=np.complex64)
#phase[floor(diffSet.shape[1]/2 - probeInnerSize/2) :floor(diffSet.shape[1]/2 - probeInnerSize/2)+probeInnerSize,  floor(diffSet.shape[2]/2 - probeInnerSize/2):floor(diffSet.shape[2]/2 - probeInnerSize/2)+probeInnerSize] = np.pi/2

## Create square phase 
#phase = np.pi * np.ones((diffSet.shape[1],diffSet.shape[2]), dtype=np.complex64)
# create complex probe
probe = probe * np.exp(1j*phase)

# size of one pixel in objectplane. (blir annorlunda för att Nx och Ny är olika)
xpixel = xfactor/pixel_det
ypixel = yfactor/pixel_det

# what the width of the diffraction pattern equals to in object plan (pixlar * pixelstorlekx/y)
sizeDiffObjectx =  Nx * xpixel
sizeDiffObjecty =  Ny * ypixel

# hur långt motorn rör sig i x och yled: 
motorWidthx = ( motorpositionx.max() - motorpositionx.min() ) * 1E-6
motorWidthy = ( motorpositiony.max() - motorpositiony.min() ) * 1E-6

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
positiony = (motorpositiony - motorpositiony.min() ) *1E-6
positionx = (motorpositionx - motorpositionx.min() ) *1E-6

## ska dessa användas ? vet ej men får samma (!?) resultat
#positiony = abs(motorpositiony - motorpositiony.max() ) *1E-6
#positionx = abs(motorpositionx - motorpositionx.max() ) *1E-6
#
## mirror diffraction patterns
#diffSet = np.fliplr(diffSet)
#
#plt.figure()                                #                        x                  y
#plt.imshow(abs(probe), cmap='gray', interpolation='none', extent=[0,sizeDiffObjecty*1E6,0,sizeDiffObjecty*1E6])
#plt.xlabel(' [µm]')
#plt.ylabel(' [µm]')
#plt.title('Initial probe amplitude')
#plt.colorbar()
#
#plt.figure()
#plt.imshow(np.angle(probe), cmap='gray', interpolation='none', extent=[0,sizeDiffObjecty*1E6,0,sizeDiffObjecty*1E6])
#plt.xlabel(' [µm]')
#plt.ylabel(' [µm]')
#plt.title('Initial probe phase')
#plt.colorbar()

# run ePIE
objectFunc, probe, ani, sse, psi = ePIE(diffSet, probe, objectFuncNy, objectFuncNx, ypixel, xpixel, positiony, positionx, nbr_scans)

#TODO someting.. PRTF?
### make ePIE function return psi (exit wave) for every position of probe
# sum over all positions?
#psi = sum(psi)
#fft1 = fft.fft(psi)
#plt.figure()
#plt.plot(fft1)

# function creates a gaussian  with amplitude A, center of function c, and sigma
def gauss(x, A, c, sigma):
    return A*np.exp(-(x-c)**2/(2*sigma**2))

xCol = np.linspace(0,probe.shape[1]-1,probe.shape[1])
xRow = np.linspace(0,probe.shape[0]-1,probe.shape[0])     #fler punkter?
yFit = gauss(xCol,152,95,2)     #sumcolumns

yCol_data = abs(probe.sum(axis=0))
yRow_data = abs(probe.sum(axis=1))
# fit probe columns to gaussian #p0 initial guesses for fitting (optional)
poptCol, puCol = curve_fit(gauss, xCol, yCol_data, p0=[yCol_data.max(),95,1])
# fit probe rows to gaussian
poptRow, puRow = curve_fit(gauss, xRow, yRow_data )

FWHM_col = 4.29193 * poptCol[2]      #(=2 * (sqrt(2*ln(10) ) )))
FWHM_row = 4.29193 * poptRow[2]
np.disp(FWHM_col)
np.disp(FWHM_row)

#2d gaussian fit
def gauss2d(xytuple, A, cx, cy, sigmax, sigmay):
    (x,y) = xytuple    # hur funkar detta?
    g = A*np.exp(- ((x-cx)**2 /(2*sigmax**2) + (y-cy)**2 /(2*sigmay**2)   ))
    return g.ravel()

# Create x and y indices
x = np.linspace(0,probe.shape[1]-1,probe.shape[1])
y = np.linspace(0,probe.shape[0]-1,probe.shape[0])
x, y = np.meshgrid(x, y)
xytuple = (x,y);
y2sgauss = gauss2d(xytuple, abs(probe).max(), 82, 95, 1, 1 )

data2d = abs(probe)
popt2d, pu2 = curve_fit(gauss2d, xytuple, data2d.ravel(), p0=[data2d.max(), 92, 81, 1, 1] )
#
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
plt.show() #show animation


# Make a centred line in x and y direction on the diffraction patterns
#diffSet[:,:,int(diffSet.shape[2]/2)] = 1 
#diffSet[:,int(diffSet.shape[1]/2),:] = 1

###plot the trimmed and centered sum of diffPatterns
#linex = np.linspace(0,diffSet.shape[2],diffSet.shape[2]+1)
#liney = np.linspace(0,diffSet.shape[1],diffSet.shape[1]+1)
#lineyy = diffSet.shape[1]/2 * np.ones((linex.shape))
def plot():
    plt.figure()
    plt.imshow(np.log10(sum(diffSet)))
    #plt.imshow(nollor)
    #plt.plot( lineyy, linex)
    plt.title('All diffraction patterns summed')
    plt.colorbar()
    
    
    #def plott():  
    plt.figure()     #, origin="lower"                         # sets the scale on axes. Should be calculated for every new experiment
    plt.imshow( np.angle(objectFunc), interpolation='none', extent=[0,objectFuncNy*ypixel*1E6, 0,objectFuncNx*xpixel*1E6])
    #plt.gca().invert_yaxis() 
    plt.xlabel(' [µm]')
    plt.ylabel(' [µm]')
    plt.title('Object phase')
    #plt.clim(-np.pi,np.pi)
    plt.colorbar()
       
    plt.figure()                                                            # horisontalt vertikalt. xpixel * size(objectfunc[xled])
    plt.imshow(abs(objectFunc), cmap='gray', interpolation='none', extent=[0,objectFuncNy*ypixel*1E6, 0,objectFuncNx*xpixel*1E6])
    plt.xlabel(' [µm]')
    plt.ylabel(' [µm]')
    plt.title('Object amplitude')
    plt.colorbar()
        
    plt.figure()
    plt.imshow(abs(probe), interpolation='none', extent=[0,sizeDiffObjecty*1E6, 0,sizeDiffObjectx*1E6])
    plt.xlabel(' [µm]')
    plt.ylabel(' [µm]')
    plt.title('Probe amplitude')
    plt.colorbar()
    
    plt.figure()                                                            # horisontalt vertikalt
    plt.imshow(np.angle(probe), interpolation='none', extent=[0,sizeDiffObjecty*1E6, 0,sizeDiffObjectx*1E6])
    plt.xlabel(' [µm]')
    plt.ylabel(' [µm]')
    plt.title('Probe phase')
    #plt.clim(-np.pi,np.pi)
    plt.colorbar()
    
    plt.figure()
    plt.plot(sse)
    plt.xlabel(' iterations ')
    plt.ylabel(' SSE ')
    plt.title('SSE')
    
    plt.figure()
    plt.plot(abs(probe.sum(axis=0)), 'b+:', label='data')                                    
    plt.plot(xCol, gauss(xCol, *poptCol), 'r-', label='fit')
    plt.plot(xCol, yFit, 'g', label='manual fit')
    plt.xlabel(' [µm]')
    plt.title('Probe summed over all columns')
    plt.legend()
    #plt.axis((0,xpixel*diffSet.shape[2]*1E6,0,10))
    
    #my_xticks = xpixel*diffSet.shape[2]*1E6
    #plt.set_xticks( my_xticks )
    
    plt.figure()
    plt.plot(abs(probe.sum(axis=1)), 'b+:', label='data')  #sum
    plt.plot(xRow,gauss(xRow, *poptRow),'r-', label='fit')
    plt.legend() 
    #plt.axis([0, sizeDiffObjectx*1E6, 0, sizeDiffObjecty*1E6])    #kontrollera så att x är x och y är y 
    #plt.yscale( )
    #plt.axis.set_xscale(sizeDiffObjectx*1E6) 
    plt.xlabel(' [µm]')
    plt.title('Probe summed over all rows')

    return 0
#plot()