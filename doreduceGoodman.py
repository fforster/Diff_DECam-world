#!/usr/bin/python2.7

import os
import re # use regular patterns
import sys, getopt # system commands
import string # string functions4
import math
import numpy as np # numerical tools
from scipy import linalg as scipylinalg
from scipy import stats
from scipy import ndimage
from scipy import signal as scipysignal
from pylab import *
from mx.DateTime import * # date conversion

from projection import projfast
from conv2fast import conv2fast

import minuit

import pyfits as fits

projfast.set_num_threads(4)
conv2fast.set_num_threads(4)

download = False
dodetrend = False
domosaic = False
doconvolve = False
dophotometry = False
obsdate = '20140324'

try:
    opts, args = getopt.getopt(sys.argv[1:], 'hfsondtmcp', ['help', 'filter=', 'supernova=', 'obsdate=', 'order=', 'download', 'detrend', 'mosaic', 'convolve', 'photometry'])
except getopt.GetoptError:
    print 'doreduce.py --help'
    sys.exit(1)
for opt, arg in opts:
    if opt in ('-h', '--help'):
        print "\nCC by F.Forster\nMatch DECam and SOI images to compute image subtraction"
        print "-f --filter: select filter"
        print "-s --supernova: select supernova"
        print "-o --obsdate: observation date (directory name)"
        print "-n --order: order of coordinate transformation (1 or 2 allowed)"
        print "-d --download: download required images according looking for MJD"
        print "-t --detrend: do bias and flat fielding"
        print "-m --mosaic: create SOI image mosaic based on the available data"
        print "-c --convolve: compute kernel, projection and image convolution"
        print "-p --photometry: do photometry"
        print "e.g. python doreduceSOI.py --filter g --supernova Bel --obsdate 20140324 --order 3 --download --detrend --mosaic --convolve --photometry"
        sys.exit()
    elif opt in ('-f', '--filter'):
        filter = arg
    elif opt in ('-s', '--supernova'):
        supernova = arg
    elif opt in ('-o', '--obsdate'):
        obsdate = arg
    elif opt in ('-n', '--order'):
        order = int(arg)
        if not order in (0, 1, 2, 3):
            print "WARNING: coordinate transformation order must be 1, 2 or 3"
            sys.exit()
    elif opt in ('-d', '--download'):
        download = True
    elif opt in ('-t', '--detrend'):
        dodetrend = True
    elif opt in ('-m', '--mosaic'):
        domosaic = True
    elif opt in ('-c', '--convolve'):
        doconvolve = True
    elif opt in ('-p', '--photometry'):
        dophotometry = True

print supernova

# other options
dotrim = True
DECamize = True  # rotate image to have the same orientation of the DECam images
docrblaster = True
doplot = True

# background backfilter size
backsize = 128

# small tick labels
matplotlib.rc('xtick', labelsize = 7) 
matplotlib.rc('ytick', labelsize = 7) 

# DECam pixel scale in SOAR pixel units
npix = 0.265 / 0.155

# open and plot reference image

refdir = '/home/fforster/Work/DATA_CMMPIPE'
refdir = '/media/fforster_fondecyt_2Tb_2014/DATA/DATA_CMMPIPE'
if supernova == 'Bel':
    field = 'Blind14A-P_10'
    CCD = 'N6'
    epoch = '02'
    coords = [2808, 986]
elif supernova == 'Dana':
    field = 'Blind14A-P_08'
    CCD = 'N28'
    epoch = '02'
    coords = [274, 352]
elif supernova == 'Greta':
    field = 'Blind14A-P_30'
    CCD = 'S31'
    epoch = '02'
    coords = [211, 1308]
elif supernova == 'Elise':
    field = 'Blind14A-P_04'
    CCD = 'S23'
    epoch = '02'
    coords = [1682, 1588]
elif supernova == 'Emilia':
    field = 'Blind14A-P_37'
    CCD = 'N30'
    epoch = '02'
    coords = [1187, 1705]
elif supernova == 'Kora':
    field = 'Blind14A-P_01'
    CCD = 'S23'
    epoch = '03'
    coords = [1840, 2159]
elif supernova == 'Laura':
    field = 'Blind14A-P_27'
    CCD = 'N31'
    epoch = '02'
    coords = [107, 1466]
elif supernova == 'Mara':
    field = 'Blind14A-P_07'
    CCD = 'N17'
    epoch = '02'
    coords = [848, 721]
elif supernova == 'Pamela':
    field = 'Blind14A-P_16'
    CCD = 'S13'
    epoch = '02'
    coords = [3793, 597]
elif supernova == 'Catalina':
    field = 'Blind14A-P_01'
    CCD = 'S19'
    epoch = '03'
    coords = [2726, 1574]
elif supernova == 'Farah':
    field = 'Blind14A-P_26'
    CCD = 'N11'
    epoch = '02'
    coords = [2744, 388]
elif supernova == 'Milena':
    field = 'Blind14A-P_09'
    CCD = 'N4'
    epoch = '02'
    coords = [1021, 201]
elif supernova == 'Nala':
    field = 'Blind14A-P_08'
    CCD = 'N21'
    epoch = '02'
    coords = [1407, 113]
else:
    print "WARNING: supernova %s not in the list of valid supernova" % supernova
    sys.exit()

# SOAR CCD which we want to project
CCDSOAR = 2

# filters
filters = {'g': 's0011 g SDSS', 'r': 's0012 r SDSS'}

# output directory
indir = '/media/fforster_fondecyt_2Tb_2014/DATA/SOI/%s' % obsdate
outdir = "%s/OUT" % indir
if not os.path.exists(outdir):
    os.makedirs(outdir)

# download DECam data
if download:

    print "Downloading DECam data to local disk..."
    command = "python ./filldata.py %s %s %s" % (field, CCD, epoch)
    print command
    os.system(command)

# fits files
files = os.listdir(indir)

# file calibration
if dodetrend:

    # initialize bias and flats
    for i in range(4):
        exec("zero%i = None" % (i + 1))
        exec("flat%i = None" % (i + 1))

    # compute bias

    print "\nComputing master bias...\n"

    for i in sorted(files):
    
        filei = "%s/%s" % (indir, i)
        if re.match("zero.*?", i):
            print "Bias", i
            obs = fits.open(filei)
            if zero1 == None:
                for iCCD in range(4):
                    zero = obs[iCCD + 1].data
                    if dotrim:
                        datasec = obs[iCCD + 1].header['DATASEC']
                        (l1, l2, l3, l4) = re.findall("\[(\d+):(\d+),(\d+):(\d+)\]", datasec)[0]
                        zero = zero[int(l3) - 1: int(l4), int(l1) - 1: int(l2)]
                    if DECamize:
                        zero = zero.transpose()
                        zero = zero[::-1,::-1]
                    exec("zero%i = zero" % (iCCD + 1))
            else:
                for iCCD in range(4):
                    zero = obs[iCCD + 1].data
                    if dotrim:
                        datasec = obs[iCCD + 1].header['DATASEC']
                        (l1, l2, l3, l4) = re.findall("\[(\d+):(\d+),(\d+):(\d+)\]", datasec)[0]
                        zero = zero[int(l3) - 1: int(l4), int(l1) - 1: int(l2)]
                    if DECamize:
                        zero = zero.transpose()
                        zero = zero[::-1,::-1]
                    exec("zero%i = np.dstack([zero%i, zero])" % (iCCD + 1, iCCD + 1))
    
    # take median bias and save it
    for iCCD in range(4):
    
        exec("zero%i = np.average(zero%i, axis = 2)" % (iCCD + 1, iCCD + 1))
        exec("obs[iCCD + 1].data = zero%i" % (iCCD + 1))
    
        exec("header%i = obs[iCCD + 1].header" % (iCCD + 1))
        exec("obs%i = fits.PrimaryHDU(data = float32(zero%i), header = header%i)" % (iCCD + 1, iCCD + 1, iCCD + 1))
        exec("obs%i.writeto(\"%s/bias_master_%02i.fits\", clobber = True)" % (iCCD + 1, outdir, iCCD + 1))
    
    # compute flats

    print "\nComputing master flats...\n"

    for i in sorted(files):

        filei = "%s/%s" % (indir, i)
        if re.match("sflat.*?", i) or re.match("flat.*?", i):
            filteri1 = fits.open(filei)[0].header['FILTER1']
            filteri2 = fits.open(filei)[0].header['FILTER2']
            if filteri1 == "s0000 Open" and filteri2 == filters[filter]:
                print "Flat", i, filter
                obs = fits.open(filei)
                if flat1 == None:
                    for iCCD in range(4):
                        flat = obs[iCCD + 1].data
                        # remove completely saturated frames
                        if np.median(flat) == np.max(flat):
                            print "Skipping saturated flat..."
                            continue
                        exptime = float(obs[0].header['EXPTIME'])
                        if dotrim:
                            datasec = obs[iCCD + 1].header['DATASEC']
                            (l1, l2, l3, l4) = re.findall("\[(\d+):(\d+),(\d+):(\d+)\]", datasec)[0]
                            flat = flat[int(l3) - 1: int(l4), int(l1) - 1: int(l2)]
                        if DECamize:
                            flat = flat.transpose()
                            flat = flat[::-1,::-1]
                        exec("flat%i = (flat - zero%i) / exptime" % (iCCD + 1, iCCD + 1))
                else:
                    for iCCD in range(4):
                        flat = obs[iCCD + 1].data
                        # remove completely saturated frames
                        if np.median(flat) == np.max(flat):
                            print "Skipping saturated flat..."
                            continue
                        if dotrim:
                            datasec = obs[iCCD + 1].header['DATASEC']
                            (l1, l2, l3, l4) = re.findall("\[(\d+):(\d+),(\d+):(\d+)\]", datasec)[0]
                            flat = flat[int(l3) - 1: int(l4), int(l1) - 1: int(l2)]
                        if DECamize:
                            flat = flat.transpose()
                            flat = flat[::-1,::-1]
                        exec("flat%i = np.dstack([flat%i, flat])" % (iCCD + 1, iCCD + 1))
    
    # take median flat and save
    for iCCD in range(4):
        exec("flat%i = np.sum(flat%i, axis = 2)" % (iCCD + 1, iCCD + 1))
        exec("flat%i = flat%i / np.median(flat%i)" % (iCCD + 1, iCCD + 1, iCCD + 1))
        exec("obs[iCCD + 1].data = flat%i" % (iCCD + 1))
    
        exec("header%i = obs[iCCD + 1].header" % (iCCD + 1))
        exec("obs%i = fits.PrimaryHDU(data = float32(flat%i), header = header%i)" % (iCCD + 1, iCCD + 1, iCCD + 1))
        exec("obs%i.writeto(\"%s/flat_master_%s_%02i.fits\", clobber = True)" % (iCCD + 1, outdir, filter, iCCD + 1))
    
    # correct images                    

    print "\nApplying bias and flat field corrections to raw images...\n"

    for i in sorted(files):

        filei = "%s/%s" % (indir, i)
        if re.match("image.*?", i):
    
            ifile = int(i[6:9])

            obs = fits.open(filei)
            header = obs[0].header
    
            filteri1 = header['FILTER1']
            filteri2 = header['FILTER2']
            obj = header['OBJECT']

            #print i, filteri1, filteri2, obj, filteri1 == "s0000 Open", filteri2 == filters[filter]
    
            if filteri1 == "s0000 Open" and filteri2 == filters[filter]:
     
                if obj != supernova:
                    continue
    
                print "\nImage", i, filter, "obj:", obj

                for iCCD in range(4):
                    
                    if iCCD + 1 != CCDSOAR:
                        continue

                    image = obs[iCCD + 1].data
                    if dotrim:
                        datasec = obs[iCCD + 1].header['DATASEC']
                        (l1, l2, l3, l4) = re.findall("\[(\d+):(\d+),(\d+):(\d+)\]", datasec)[0]
                        image = image[int(l3) - 1: int(l4), int(l1) - 1: int(l2)]
                    if DECamize:
                        image = image.transpose()
                        image = image[::-1,::-1]
                    exec("image%i = (image - zero%i) / flat%i" % (iCCD + 1, iCCD + 1, iCCD + 1))
                    exec("obs[iCCD + 1].data = image%i" % (iCCD + 1))
                    exec("header = obs[%i].header" % (iCCD + 1))
                    header.update('GAIN', 2.0)
                    delfields = ['NEXTEND', 'IMAGEID', 'DETSIZE', 'OFTDIR', 'WAT0_001', 'DATASEC', 'BIASSEC', 'CCDSEC', 'AMPSEC', 'TRIMSEC', 'DETSEC', 'CCDSIZE', 'NCCDS', 'NAMPS', 'TV1FOC', 'IMAGESWV', 'OFTDIR', 'XTALKFIL', 'RECNO', 'PREFLASH', 'RADECEQ', 'CCDNAME', 'AMPNAME', 'RDNOISE', 'SATURATE', 'AMPINTEG', 'BPM', 'WCSASTRM', 'EQUINOX', 'TELEQIN', 'WAT1_001', 'WAT1_002', 'WAT1_003', 'WAT1_004', 'WAT1_005', 'WAT2_001', 'WAT2_002', 'WAT2_003', 'WAT2_004', 'WAT2_005', 'CHECKSUM', 'DATASUM', 'CHECKVER']
                    for delfield in delfields:
                        try:
                            del header[delfield]
                        except:
                            print "Cannot delete field %s" % delfield
                    header.update('CTYPE1', 'RA--TAN')
                    header.update('CTYPE2', 'DEC--TAN')
                    exec("obs%i = fits.PrimaryHDU(data = float32(image%i), header = header)" % (iCCD + 1, iCCD + 1))
                    filename = "%s/%s_%s_%02i_%04i.fits" % (outdir, obj, filter, iCCD + 1, ifile)
                    exec("obs%i.writeto(\"%s\", clobber = True)" % (iCCD + 1, filename))

                    # run sextractor
                    background = filename.replace(".fits", "_background-%03i.fits" % backsize)
                    command = "sex %s -CATALOG_NAME %s-catalogue.dat -BACK_SIZE %i -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME %s -VERBOSE_TYPE QUIET" % (filename, filename, backsize, background)
                    print command
                    os.system(command)
    
                    # load catalogue
                    try:
                        (x, y, flux, e_flux, r, flag) = np.loadtxt("%s-catalogue.dat" % (filename), usecols = (1, 2, 5, 6, 8, 9)).transpose()
                    except:
                        print "Sextractor failed"
                        continue
    
                    # plot image + catalogue
                    if doplot:
                        fig, ax = plt.subplots(3, 1, figsize = (17, 15))
                        fig.subplots_adjust(wspace = 0, hspace = 0)
                        
                        ax[0].scatter(x[flag <= 2], y[flag <= 2], marker = 'o', facecolors = 'none', color = 'white', s = 20 * r)
                        ax[1].scatter(x[flag <= 2], y[flag <= 2], marker = 'o', facecolors = 'none', color = 'white',  s = 20 * r)
                        ax[2].scatter(x[flag <= 2], y[flag <= 2], marker = 'o', facecolors = 'none', color = 'white', s = 20 * r)
                        
                        ax[0].scatter(x[flag <= 4], y[flag <= 4], marker = 'o', facecolors = 'none', color = 'y', s = 20 * r)
                        ax[1].scatter(x[flag <= 4], y[flag <= 4], marker = 'o', facecolors = 'none', color = 'y', s = 20 * r)
                        ax[2].scatter(x[flag <= 4], y[flag <= 4], marker = 'o', facecolors = 'none', color = 'y', s = 20 * r)
                        
                        ax[1].scatter(x[flag > 4], y[flag > 4], marker = 'o', facecolors = 'none', color = 'r', s = 20 * r)
                        ax[2].scatter(x[flag > 4], y[flag > 4], marker = 'o', facecolors = 'none', color = 'r', s = 20 * r)
                        ax[0].scatter(x[flag > 4], y[flag > 4], marker = 'o', facecolors = 'none', color = 'r', s = 20 * r)
                        
                        image = fits.open(filename)[0].data
                        ax[0].imshow(image, cmap = 'gray', interpolation = 'nearest', clim = (np.percentile(image.flatten(), 10), np.percentile(image.flatten(), 99)), origin = 'lower')
                        bg = fits.open(background)[0].data
                        image = image - bg
                        ax[1].imshow(image, cmap = 'gray', interpolation = 'nearest', clim = (np.percentile(image.flatten(), 10), np.percentile(image.flatten(), 99)), origin = 'lower')
                        ax[2].imshow(bg, cmap = 'gray', interpolation = 'nearest', origin = 'lower')
                        plt.savefig("%s/%s_%s_%02i_%04i_nosky-%03i_o%i.png" % (outdir, obj, filter, iCCD + 1, ifile, backsize, order), bbox_inches = 'tight', dpi = 150)

else:
    
    "Skipping Image detrend.."


# PROJECTION
# ----------------------------------------

# find linear transformation relating two sets of points
def findshifttransformation(x1, y1, x2, y2):

    # solve linear transformation between two coordinate systems
    # find best transformation relating all these points
    # need to write the system of equations x' = ax + b as X' = X beta
    # we use beta = (a11 b1 a22 b2)^T
    # where x' = a11 x + b1 and y' = a22 y + b2 for npt points
    # then X' = (x1...xn y1...yn)^T, X = ((x1 1 0 0) (x2 1 0 0) ... (xn 1 0 0) (0 0 y1 1) (0 0 y2 1) ... (0 0 yn 1))
    # and (X^T X) beta = (X^T X')

    npt = len(x1)
    if len(x1new) < 2:
        print "\n\nWARNING: Not enough stars to shift astrometric solution (%i)...\n\n" % (len(x1new))
        return None
    Y = np.zeros(2 * npt)
    Y[0:npt] = x2
    Y[npt: 2 * npt] = y2
    X = np.zeros((2. * npt, 4))
    X[0:npt, 0] = x1
    X[0:npt, 1] = 1.
    X[npt: 2 * npt, 2] = y1
    X[npt: 2 * npt, 3] = 1.
    # solve
    mat = np.dot(X.transpose(), X)
    rhs = np.dot(X.transpose(), Y)
    try:
        print "Solving shift+scale transformation..."
        (a11, b1, a22, b2) = scipylinalg.solve(mat, rhs)
    except:
        print "Cannot solve linear system"
        return None

    return (a11, b1, a22, b2)

# find linear transformation relating two sets of points
def findlineartransformation(x1, y1, x2, y2):

    # solve linear transformation between two coordinate systems
    # find best transformation relating all these points
    # need to write the system of equations x' = ax + b as X' = X beta
    # we use beta = (a11 a12 b1 a21 a22 b2)^T
    # where x' = a11 x + a12 y + b1 and y' = a21 x + a22 y + b2 for npt points
    # then X' = (x1...xn y1...yn)^T, X = ((x1 y1 1 0 0 0) (x2 y2 1 0 0 0) ... (xn yn 1 0 0 0) (0 0 0 x1 y1 1) (0 0 0 x2 y2 1) ... (0 0 0 xn yn 1))
    # and (X^T X) beta = (X^T X')

    npt = len(x1)
    if len(x1new) < 3:
        print "\n\nWARNING: Not enough stars to do linear astrometric solution (%i)...\n\n" % (len(x1new))
        return None
    Y = np.zeros(2 * npt)
    Y[0:npt] = x2
    Y[npt: 2 * npt] = y2
    X = np.zeros((2. * npt, 6))
    X[0:npt, 0] = x1
    X[0:npt, 1] = y1
    X[0:npt, 2] = 1.
    X[npt: 2 * npt, 3] = x1
    X[npt: 2 * npt, 4] = y1
    X[npt: 2 * npt, 5] = 1.
    # solve
    mat = np.dot(X.transpose(), X)
    rhs = np.dot(X.transpose(), Y)
    try:
        print "Solving linear transformation..."
        (a11, a12, b1, a21, a22, b2) = scipylinalg.solve(mat, rhs)
    except:
        print "Cannot solve linear system"
        return None

    return (a11, a12, b1, a21, a22, b2)

# find non linear transformation relating two sets of points
def find2ndordertransformation(x1, y1, x2, y2):

    # solve 2nd order transformation between two coordinate systems
    # find best transformation relating all these points
    # need to write the system of equations:
    # x' = a1 x + b11 x + b12 y + c1 xy + d11 x^2 + d12 y^2
    # y' = a2 y + b21 x + b22 y + c2 xy + d21 x^2 + d22 y^2
    # X' = X beta
    # we use beta = (a1 b11 b12 c1 d11 d12 a2 b21 b22 c2 d21 d22)^T
    # then X' = (x1...xn y1...yn)^T, X = ((1 x1 y1 x1*y1 x1^2 y1^2 0 0 0 0 0 0) ... (1 x1 y1 x1*y1 x1^2 y1^2 0 0 0 0 0 0) (0 0 0 0 0 0 1 x1 y1 x1*y1 x1^2 y1^2) ... (0 0 0 0 0 0 1 xn yn xn*yn xn^2 yn^2)
    # the least squares errors is found that beta which is solution of the following linear system
    # (X^T X) beta = (X^T X')
    # below we use the notation X'->Y

    npt = len(x1)
    if len(x1new) < 6:
        print "\n\nWARNING: Not enough stars to do quadratic astrometric solution (%i)...\n\n" % (len(x1new))
        return None
    Y = np.zeros(2 * npt)
    Y[0:npt] = x2
    Y[npt: 2 * npt] = y2
    X = np.zeros((2. * npt, 12))
    X[0: npt, 0] = 1.
    X[0: npt, 1] = x1
    X[0: npt, 2] = y1
    X[0: npt, 3] = x1 * y1
    X[0: npt, 4] = x1 * x1
    X[0: npt, 5] = y1 * y1
    X[npt: 2 * npt, 6] = X[0:npt, 0]
    X[npt: 2 * npt, 7] = X[0:npt, 1]
    X[npt: 2 * npt, 8] = X[0:npt, 2]
    X[npt: 2 * npt, 9] = X[0:npt, 3]
    X[npt: 2 * npt, 10] = X[0:npt, 4]
    X[npt: 2 * npt, 11] = X[0:npt, 5]
    # solve
    mat = np.dot(X.transpose(), X)
    rhs = np.dot(X.transpose(), Y)
    try:
        print "Solving 2nd order transformation..."
        (a1, b11, b12, c1, d11, d12, a2, b21, b22, c2, d21, d22) = scipylinalg.solve(mat, rhs)
    except:
        print "Cannot solve linear system"
        return None

    return (a1, b11, b12, c1, d11, d12, a2, b21, b22, c2, d21, d22)


# find non linear transformation relating two sets of points
def find3rdordertransformation(x1, y1, x2, y2):

    # solve 3rd order transformation between two coordinate systems
    # find best transformation relating all these points
    # need to write the system of equations:
    # x' = a1 x + b11 x + b12 y + c1 xy + d11 x^2 + d12 y^2 + e1 * x^2y + f1 * xy^2 + g11 * x^3 + g12 * y^3
    # y' = a2 y + b21 x + b22 y + c2 xy + d21 x^2 + d22 y^2 + e2 * x^2y + f2 * xy^2 + g21 * x^3 + g22 * y^3
    # X' = X beta
    # we use beta = (a1 b11 b12 c1 d11 d12 e1 f1 g11 g12 a2 b21 b22 c2 d21 d22 e2 f2 g21 g22)^T
    # then X' = (x1...xn y1...yn)^T, X = ((1 x1 y1 x1*y1 x1^2 y1^2 x1^2*y1 x1*y1^2 x1^3 y1^3 0 0 0 0 0 0 0 0 0 0) ... (1 xn yn xn*yn xn^2 yn^2 xn^2*yn xn*yn^2 xn^3 yn^3 0 0 0 0 0 0 0 0 0) (0 0 0 0 0 0 0 0 0 0 1 x1 y1 x1*y1 x1^2 y1^2 x1^2*y1 x1*y1^2 x1^3 y1^3) ... (0 0 0 0 0 0 0 0 0 0 1 xn yn xn*yn xn^2 yn^2 xn^2*yn xn*yn^2 xn^3 yn^3)
    # the least squares errors is found that beta which is solution of the following linear system
    # (X^T X) beta = (X^T X')
    # below we use the notation X'->Y

    npt = len(x1)
    if len(x1new) < 10:
        print "\n\nWARNING: Not enough stars to do cubic astrometric solution (%i)...\n\n" % (len(x1new))
        return None
    Y = np.zeros(2 * npt)
    Y[0:npt] = x2
    Y[npt: 2 * npt] = y2
    X = np.zeros((2. * npt, 20))
    X[0: npt, 0] = 1.
    X[0: npt, 1] = x1
    X[0: npt, 2] = y1
    X[0: npt, 3] = x1 * y1
    X[0: npt, 4] = x1 * x1
    X[0: npt, 5] = y1 * y1
    X[0: npt, 6] = x1 * x1 * y1
    X[0: npt, 7] = x1 * y1 * y1
    X[0: npt, 8] = x1 * x1 * x1
    X[0: npt, 9] = y1 * y1 * y1
    X[npt: 2 * npt, 10] = X[0:npt, 0]
    X[npt: 2 * npt, 11] = X[0:npt, 1]
    X[npt: 2 * npt, 12] = X[0:npt, 2]
    X[npt: 2 * npt, 13] = X[0:npt, 3]
    X[npt: 2 * npt, 14] = X[0:npt, 4]
    X[npt: 2 * npt, 15] = X[0:npt, 5]
    X[npt: 2 * npt, 16] = X[0:npt, 6]
    X[npt: 2 * npt, 17] = X[0:npt, 7]
    X[npt: 2 * npt, 18] = X[0:npt, 8]
    X[npt: 2 * npt, 19] = X[0:npt, 9]
    # solve
    mat = np.dot(X.transpose(), X)
    rhs = np.dot(X.transpose(), Y)
    try:
        print "Solving 3rd order transformation (npt: %i)..." % npt
        (a1, b11, b12, c1, d11, d12, e1, f1, g11, g12, a2, b21, b22, c2, d21, d22, e2, f2, g21, g22) = scipylinalg.solve(mat, rhs)
    except:
        print "Cannot solve linear system"
        return None

    return (a1, b11, b12, c1, d11, d12, e1, f1, g11, g12, a2, b21, b22, c2, d21, d22, e2, f2, g21, g22)


# return intersection of two samples given deltas
def intersect(x1, x2, y1, y2, deltax, deltay, npix):
    
    x1new = []
    x2new = []
    y1new = []
    y2new = []

    for isource in range(len(x2)):
        
        distx = (x2[isource] / npix + deltax - x1)
        disty = (y2[isource] / npix + deltay - y1)
        dist = distx * distx + disty * disty
        idx = np.argmin(dist)

        maxdist = 5
        if dist[idx] < maxdist**2:
            x1new.append(x1[idx])
            x2new.append(x2[isource])
            y1new.append(y1[idx])
            y2new.append(y2[isource])
    
    x1new = np.array(x1new)
    x2new = np.array(x2new)
    y1new = np.array(y1new)
    y2new = np.array(y2new)

    return x1new, y1new, x2new, y2new
    
# distance between two sets
def xydistance(x1, x2, y1, y2, delta):
    
    nsources = 0
    total = 0

    for isource in range(len(x2)):
        
        distx = (x2[isource] - x1)
        disty = (y2[isource] - y1)
        dist = distx * distx + disty * disty
        idx = np.argmin(dist)

        if dist[idx] < delta**2:
            nsources += 1
            total += dist[idx]
    
    return nsources, total

# function that tries to find first course astrometric solution
def roughastro(x1, x2, y1, y2, deltaxmin, deltaxmax, deltaymin, deltaymax, delta):

    ibest = 0
    jbest = 0
    nbest = 0

    for i in arange(deltaxmin, deltaxmax, delta):
        
        for j in arange(deltaymin, deltaymax, delta):

            (nsources, dist) = xydistance(x1, x2 + i, y1, y2 + j, delta)
            if nsources >= nbest:
                ibest = i
                jbest = j
                nbest = nsources
                #print ibest, jbest, nsources, dist
    
    return ibest, jbest

# open reference file
fitsref = "%s/%s/%s/%s_%s_%s_image_crblaster.fits" % (refdir, field, CCD, field, CCD, epoch)
headerref = fits.open(fitsref)[0].header
dataref = fits.open(fitsref)[0].data

MJDref = headerref['MJD-OBS']

# load reference catalogue
(xref, yref, fluxref, e_fluxref, rref, flagref) = np.loadtxt("%s-catalogue_wtmap.dat" % (fitsref), usecols = (1, 2, 5, 6, 8, 9)).transpose()

# compute array sizes
(nx1, ny1) = np.shape(dataref)

if domosaic:

    final = np.zeros(np.shape(dataref))
    finalbg = np.zeros(np.shape(dataref))
    nmosaic = np.zeros(np.shape(dataref), dtype = int)
    
    ifiles = []
    iMJDs = []

    for i in sorted(files):
            
        if re.match("image.*?", i):
            
            filei = "%s/%s" % (indir, i)
            ifile = int(i[6:9])
    
            obs = fits.open(filei)
            header = obs[0].header
        
            filteri1 = header['FILTER1']
            filteri2 = header['FILTER2']
            obj = header['OBJECT']
            MJD = float(header['MJD-OBS'])
    
            #print i, filteri1, filteri2, obj, filteri1 == "s0000 Open", filteri2 == filters[filter]
        
            if filteri1 == "s0000 Open" and filteri2 == filters[filter]:
         
                if obj != supernova:
                    continue
    
                print "\n\nImage", i, filter, "obj:", obj

                ifiles.append(ifile)
                iMJDs.append(MJD)

                # open image to project
                fitsnew = "%s/%s_%s_%02i_%04i.fits" % (outdir, supernova, filter, CCDSOAR, ifile)
                background = fitsnew.replace(".fits", "_background-%03i.fits" % backsize)
                headernew = fits.open(fitsnew)[0].header
                background = fits.open(background)[0].data
                datanew = fits.open(fitsnew)[0].data - background
                (nx2, ny2) = np.shape(datanew)
    
                # load sextractor sources
                (x, y, flux, e_flux, r, flag) = np.loadtxt("%s/%s_%s_%02i_%04i.fits-catalogue.dat" % (outdir, supernova, filter, CCDSOAR, ifile), usecols = (1, 2, 5, 6, 8, 9)).transpose()
    
                # plot reference image
                if doplot:
                    fig, ax = plt.subplots()
                    ax.imshow(dataref, interpolation = 'nearest', cmap = 'gray', clim = (np.percentile(dataref, 1), np.percentile(dataref, 99)), origin = 'lower')
                    ax.scatter(xref[flagref <= 2], yref[flagref <= 2], marker = 'o', facecolors = 'none', color = 'y', s = 10, lw = 0.1)
    
                # Matching
                # ---------------------
                    
                print "Looking for rough match..."
                
                # first guess
                ibest = coords[1] - ny2 / npix / 2. 
                jbest = coords[0] - nx2 / npix / 2.
    
                print "Refining solution to 25 pixels..."
                (ibest, jbest) = roughastro(xref[flagref <= 4], x[flag <= 4] / npix, yref[flagref <=4], y[flag <= 4] / npix, ibest - 300, ibest + 300, jbest - 300, jbest + 300, 25)
                
                print "Refining solution to 5 pixels..."
                (ibest, jbest) = roughastro(xref[flagref <= 4], x[flag <= 4] / npix, yref[flagref <=4], y[flag <= 4] / npix, ibest - 25, ibest + 25, jbest - 25, jbest + 25, 5)
                
                print "Refining solution to 2 pixels..."
                (ibest, jbest) = roughastro(xref[flagref <= 4], x[flag <= 4] / npix, yref[flagref <=4], y[flag <= 4] / npix, ibest - 5, ibest + 5, jbest - 5, jbest + 5, 2)
                
                if doplot:
                    ax.scatter(x[flag <= 4] / npix + ibest, y[flag <= 4] / npix + jbest, marker = 'o', facecolors = 'none', color = 'r', s = 10, lw = 0.2)
                    plt.savefig("imref.png", dpi = 400, bbox_inches = 'tight')
                    
                    
                # Linear transformation
                # ---------------------------
                    
                # find intersection
                (x1new, y1new, x2new, y2new) = intersect(xref[flagref <= 4], x[flag <= 4], yref[flagref <= 4], y[flag <= 4], ibest, jbest, npix)
                
                if doplot:
                    # plot intersection in DECam coordinates before projection
                    fig, ax = plt.subplots()
                    for ipt in range(len(x1new)):
                        ax.plot([x1new[ipt], x2new[ipt] / npix + ibest], [y1new[ipt], y2new[ipt] / npix + jbest], 'r')
                    plt.savefig("xy_match_DECam_o%i.png" % order)
    
                    # plot intersection in DECam coordinates before projection
                    fig, ax = plt.subplots()
                    ax.scatter(x1new - (x2new / npix + ibest), y1new - (y2new / npix + jbest), marker = '.') 
                    plt.savefig("xy_diff_DECam_o%i.png" % order)
    
                # solve linear transformation between two coordinate systems
                if order == 0:
                    trans = findshifttransformation(x1new, y1new, x2new, y2new)
                if order == 1:
                    trans = findlineartransformation(x1new, y1new, x2new, y2new)
                elif order == 2:
                    trans = find2ndordertransformation(x1new, y1new, x2new, y2new)
                elif order == 3:
                    trans = find3rdordertransformation(x1new, y1new, x2new, y2new)
                if trans != None:
                    if order == 0:
                        (a11, b1, a22, b2) = trans
                    if order == 1:
                        (a11, a12, b1, a21, a22, b2) = trans
                    elif order == 2:
                        (a1, b11, b12, c1, d11, d12, a2, b21, b22, c2, d21, d22) = trans
                    elif order == 3:
                        (a1, b11, b12, c1, d11, d12, e1, f1, g11, g12, a2, b21, b22, c2, d21, d22, e2, f2, g21, g22) = trans
                else:
                    continue
    
                # transformed coordinates
                if order == 0:
                    x1t = b1 + a11 * x1new
                    y1t = b2 + a22 * y1new
                elif order == 1:
                    x1t = b1 + a11 * x1new + a12 * y1new
                    y1t = b2 + a21 * x1new + a22 * y1new
                elif order == 2:
                    x1t = a1 + b11 * x1new + b12 * y1new + c1 * x1new * y1new + d11 * x1new * x1new + d12 * y1new * y1new
                    y1t = a2 + b21 * x1new + b22 * y1new + c2 * x1new * y1new + d21 * x1new * x1new + d22 * y1new * y1new
                elif order == 3:
                    x1t = a1 + b11 * x1new + b12 * y1new + c1 * x1new * y1new + d11 * x1new * x1new + d12 * y1new * y1new + e1 * x1new * x1new * y1new + f1 * x1new * y1new * y1new + g11 * x1new**3 + g12 * y1new**3
                    y1t = a2 + b21 * x1new + b22 * y1new + c2 * x1new * y1new + d21 * x1new * x1new + d22 * y1new * y1new + e2 * x1new * x1new * y1new + f2 * x1new * y1new * y1new + g21 * x1new**3 + g22 * y1new**3
                    
                # plot intersection in SOAR coordinates after projection
                if doplot:
                    fig, ax = plt.subplots()
                    for ipt in range(len(x1new)):
                        ax.plot([x1t[ipt], x2new[ipt]], [y1t[ipt], y2new[ipt]], 'r')
                    plt.savefig("xy_match_SOAR_o%i.png" % order)

                # remove outliers to refine solution
                diffx = x1t - x2new
                diffy = y1t - y2new
                mask = (diffx >= np.percentile(diffx, 5)) & (diffx <= np.percentile(diffx, 95)) & (diffy >= np.percentile(diffy, 5)) & (diffy <= np.percentile(diffy, 95))
                if order == 0 and np.sum(mask) < 2 or order == 1 and np.sum(mask) < 3 or order == 2 and np.sum(mask) < 6 or order == 3 and np.sum(mask) < 10:
                    mask = (x1t == x1t)
    
                # plot transformed DECam coordinates in SOAR coordinates
                if doplot:
                    fig, ax = plt.subplots()
                    ax.scatter(x1t - x2new, y1t - y2new, marker = '.', c = 'r')
                    ax.scatter(x1t[mask] - x2new[mask], y1t[mask] - y2new[mask], marker = 'o', c = 'b')
    
                # solve new linear transformation between two coordinate systems
                if order == 0:
                    trans = findshifttransformation(x1new[mask], y1new[mask], x2new[mask], y2new[mask])
                elif order == 1:
                    trans = findlineartransformation(x1new[mask], y1new[mask], x2new[mask], y2new[mask])
                elif order == 2:
                    trans = find2ndordertransformation(x1new[mask], y1new[mask], x2new[mask], y2new[mask])
                elif order == 3:
                    trans = find3rdordertransformation(x1new[mask], y1new[mask], x2new[mask], y2new[mask])
                if trans != None:
                    if order == 0:
                        (a11, b1, a22, b2) = trans
                    elif order == 1:
                        (a11, a12, b1, a21, a22, b2) = trans
                    elif order == 2:
                        (a1, b11, b12, c1, d11, d12, a2, b21, b22, c2, d21, d22) = trans
                    elif order == 3:
                        (a1, b11, b12, c1, d11, d12, e1, f1, g11, g12, a2, b21, b22, c2, d21, d22, e2, f2, g21, g22) = trans
                else:
                    break
    
                print "rms:", np.sqrt(np.sum((x1t[mask] - x2new[mask])**2 + (y1t[mask] - y2new[mask])**2) / len(x1t[mask]))

                # transformed coordinates
                if order == 0:
                    x1t = b1 + a11 * x1new
                    y1t = b2 + a22 * y1new
                elif order == 1:
                    x1t = b1 + a11 * x1new + a12 * y1new
                    y1t = b2 + a21 * x1new + a22 * y1new
                elif order == 2:
                    x1t = a1 + b11 * x1new + b12 * y1new + c1 * x1new * y1new + d11 * x1new * x1new + d12 * y1new * y1new
                    y1t = a2 + b21 * x1new + b22 * y1new + c2 * x1new * y1new + d21 * x1new * x1new + d22 * y1new * y1new
                elif order == 3:
                    x1t = a1 + b11 * x1new + b12 * y1new + c1 * x1new * y1new + d11 * x1new * x1new + d12 * y1new * y1new + e1 * x1new * x1new * y1new + f1 * x1new * y1new * y1new + g11 * x1new**3 + g12 * y1new**3
                    y1t = a2 + b21 * x1new + b22 * y1new + c2 * x1new * y1new + d21 * x1new * x1new + d22 * y1new * y1new + e2 * x1new * x1new * y1new + f2 * x1new * y1new * y1new + g12 * x1new**3 + g22 * y1new**3
                    
                # plot new solution
                if doplot:
                    ax.scatter(x1t[mask] - x2new[mask], y1t[mask] - y2new[mask], marker = 'o', c = 'g', s = 200, facecolors = 'none')
                    plt.savefig("xy_diff_SOAR_o%i.png" % order)
    
                # show all original points in new coordinate system
                if doplot:
                    fig, ax = plt.subplots()
                    if order == 0:
                        x1tref = b1 + a11 * xref
                        y1tref = b2 + a22 * yref
                    elif order == 1:
                        x1tref = b1 + a11 * xref + a12 * yref
                        y1tref = b2 + a21 * xref + a22 * yref
                    elif order == 2:
                        x1tref = a1 + b11 * xref + b12 * yref + c1 * xref * yref + d11 * xref * xref + d12 * yref * yref
                        y1tref = a2 + b21 * xref + b22 * yref + c2 * xref * yref + d21 * xref * xref + d22 * yref * yref
                    elif order == 3:
                        x1tref = a1 + b11 * xref + b12 * yref + c1 * xref * yref + d11 * xref * xref + d12 * yref * yref + e1 * xref * xref * yref + f1 * xref * yref * yref + g11 * xref**3 + g12 * yref**3
                        y1tref = a2 + b21 * xref + b22 * yref + c2 * xref * yref + d21 * xref * xref + d22 * yref * yref + e2 * xref * xref * yref + f2 * xref * yref * yref + g21 * xref**3 + g22 * yref**3
                    ax.scatter(x1tref, y1tref, marker = '.', c = 'r', lw = 0)
                    ax.scatter(x1t, y1t, marker = 'o', c = 'b', facecolors = 'none')
                    print np.max(xref), np.max(yref)
                    plt.savefig("xy_transformed_SOAR_o%i.png" % order)
    
                # fix nans
                datanew[np.invert(np.isfinite(datanew))] = 0
                background[np.invert(np.isfinite(background))] = 0
    
                # project image
                print "Projecting sky subtracted image..."
                alanczos = 4
                if order == 0:
                    projfast.o0_lanczos(alanczos, ny1, nx1, ny2, nx2, a11, a22, b1, b2, datanew.transpose())
                elif order == 1:
                    projfast.o1_lanczos(alanczos, ny1, nx1, ny2, nx2, a11, a12, a21, a22, b1, b2, datanew.transpose())
                elif order == 2:
                    projfast.o2_lanczos(alanczos, ny1, nx1, ny2, nx2, a1, b11, b12, c1, d11, d12, a2, b21, b22, c2, d21, d22, datanew.transpose())
                elif order == 3:
                    projfast.o3_lanczos(alanczos, ny1, nx1, ny2, nx2, a1, b11, b12, c1, d11, d12, e1, f1, g11, g12, a2, b21, b22, c2, d21, d22, e2, f2, g21, g22, datanew.transpose())
                datanewproj = np.array((projfast.imageout[0:ny1, 0:nx1]).transpose())
                print np.max(datanewproj)

                # project background
                print "Projecting background..."
                alanczos = 4
                if order == 0:
                    projfast.o0_lanczos(alanczos, ny1, nx1, ny2, nx2, a11, a22, b1, b2, background.transpose())
                elif order == 1:
                    projfast.o1_lanczos(alanczos, ny1, nx1, ny2, nx2, a11, a12, a21, a22, b1, b2, background.transpose())
                elif order == 2:
                    projfast.o2_lanczos(alanczos, ny1, nx1, ny2, nx2, a1, b11, b12, c1, d11, d12, a2, b21, b22, c2, d21, d22, background.transpose())
                elif order == 3:
                    projfast.o3_lanczos(alanczos, ny1, nx1, ny2, nx2, a1, b11, b12, c1, d11, d12, e1, f1, g11, g12, a2, b21, b22, c2, d21, d22, e2, f2, g21, g22, background.transpose())
                backgroundproj = np.array((projfast.imageout[0:ny1, 0:nx1]).transpose())
                print np.max(backgroundproj)

                # update headerref
                headerref.update('MJD-OBS', MJD)

                # save projected background
                proj = fits.PrimaryHDU(data = backgroundproj, header = headerref)
                fitsout = "%s/%s_%s_%02i_%04i_background_DECam_o%i.fits" % (outdir, supernova, filter, CCDSOAR, ifile, order)
                proj.writeto(fitsout, clobber = True)

                # save total projection
                print "Saving total projected image"
                datanewproj = np.array(backgroundproj + datanewproj)
                proj = fits.PrimaryHDU(data = datanewproj, header = headerref)
                fitsout = "%s/%s_%s_%02i_%04i_DECam_o%i.fits" % (outdir, supernova, filter, CCDSOAR, ifile, order)
                proj.writeto(fitsout, clobber = True)

                # run sextractor on projected image + background image
                if docrblaster:
                    
                    crblasterpath = '/home/fforster/Work/crblaster/crblaster'
                    
                    inmosaic = "%s/%s_%s_%02i_%04i_DECam_o%i.fits" % (outdir, supernova, filter, CCDSOAR, ifile, order)
                    outmosaic = "%s/%s_%s_%02i_%04i_DECam_crblaster_o%i.fits" % (outdir, supernova, filter, CCDSOAR, ifile, order)
                    backgroundfile = "%s/%s_%s_%02i_%04i_background_DECam_o%i.fits" % (outdir, supernova, filter, CCDSOAR, ifile, order)
                    
                    if os.path.exists(outmosaic):
                        os.system("rm -rf %s" % outmosaic)
                        
                    ncores = 4
                    command = "mpirun -n %i %s 1 1 %i %s %s" % (ncores, crblasterpath, ncores, inmosaic, outmosaic)
                    print command
                    os.system(command)
                    
                    # recover sky subtracted image
                    datanewprojcrblaster = fits.open(outmosaic)[0].data

                    # copy data excluding edges eaten out by crblaster
                    sky = np.median(datanewprojcrblaster[datanewprojcrblaster != 0])
                    print sky
                    datanewproj[datanewprojcrblaster > sky / 2.] = datanewprojcrblaster[datanewprojcrblaster > sky / 2.]

                    # save new total projection with crblaster
                    print "Saving total projected image"
                    proj = fits.PrimaryHDU(data = datanewproj, header = headerref)
                    proj.writeto(outmosaic, clobber = True)

                    # subtract sky
                    datanewproj = np.array(datanewproj - backgroundproj)

                # add to master mosaic
                final += datanewproj
                finalbg += backgroundproj
                nmosaic[datanewproj != 0] += 1

    # effective MJD
    MJD = np.average(np.array(iMJDs))
    headerref.update("MJD-OBS", "%s" % MJD)

    # save final projected image
    final[nmosaic != 0] = final[nmosaic != 0] / nmosaic[nmosaic != 0]
    final = float32(final)
    proj = fits.PrimaryHDU(data = final, header = headerref)
    fitsout = "%s/%s_%s_%02i_nosky_DECam_o%i.fits" % (outdir, supernova, filter, CCDSOAR, order)
    proj.writeto(fitsout, clobber = True)

    # save final projected background
    finalbg[nmosaic != 0] = finalbg[nmosaic != 0] / nmosaic[nmosaic != 0]
    finalbg = float32(finalbg)
    proj = fits.PrimaryHDU(data = finalbg, header = headerref)
    fitsout = "%s/%s_%s_%02i_background_DECam_o%i.fits" % (outdir, supernova, filter, CCDSOAR, order)
    proj.writeto(fitsout, clobber = True)

    # save number of images in the mosaic
    fitsout = "%s/%s_%s_%02i_nmosaic_DECam_o%i.fits" % (outdir, supernova, filter, CCDSOAR, order)
    proj = fits.PrimaryHDU(data = nmosaic, header = headerref)
    proj.writeto(fitsout, clobber = True)

# prepare for convolution
# -----------------------------------------

# load ivar array
nvar = 81
ivarf = np.loadtxt("ivar_%i.dat" % nvar, dtype = int) - 1  # 9, 57, 81, 289, 625
nf = np.shape(ivarf)[0]

# psf 
npsf = 21
npsf2 = npsf * npsf
ieq2i = np.zeros(npsf2, dtype = int)
ieq2j = np.zeros(npsf2, dtype = int)
for i in range(npsf):
    for j in range(npsf):
        ieq = i * npsf + j
        ieq2i[ieq] = i
        ieq2j[ieq] = j

# stamp sizes
dn = int((npsf + nf) / 2.)
npsfh = npsf / 2
nfh = nf / 2

# for star selection
Xstars, Ystars = np.meshgrid(np.array(range(npsf + nf)), np.array(range(npsf + nf)))
rs2Dstars = np.array(np.sqrt((Xstars - (npsf + nf - 1.) / 2.)**2 + (Ystars - (npsf + nf - 1.) / 2.)**2))

if doconvolve:

    # run sextractor on original image to remove background
    command = "sex %s -CATALOG_NAME %s-catalogue.dat -BACK_SIZE %i -CHECKIMAGE_TYPE -BACKGROUND -CHECKIMAGE_NAME %s -VERBOSE_TYPE QUIET" % (fitsref, fitsref, backsize, fitsref.replace(".fits", "_nosky.fits"))
    print command
    os.system(command)
    datareforig = np.array(dataref)
    dataref = fits.open(fitsref.replace(".fits", "_nosky.fits"))[0].data
    headerref = fits.open(fitsref.replace(".fits", "_nosky.fits"))[0].header

    # reference image variance map (must add wtmap)
    decamgain = np.sqrt(headerref['ARAWGAIN'])
    decamreadnoise = 7.
    varref = datareforig / decamgain 
    noise = fits.PrimaryHDU(data = np.sqrt(varref), header = headerref)
    noise.writeto(fitsref.replace(".fits", "_noise.fits"), clobber = True)

    # run sextractor on final mosaic image
    filename = "%s/%s_%s_%02i_nosky_DECam_o%i.fits" % (outdir, supernova, filter, CCDSOAR, order)
    headernew = fits.open(filename)[0].header
    datanew = fits.open(filename)[0].data
    backgroundnew = fits.open(filename.replace("nosky", "background"))[0].data
    nmosaic = fits.open(filename.replace("nosky", "nmosaic"))[0].data
    
    # new map noise map
    soarreadnoise = 4.4
    soargain = 2.0
    varnew = np.empty_like(datanew)
    varnew[nmosaic > 0] = (soarreadnoise**2 + (datanew[nmosaic > 0] + backgroundnew[nmosaic > 0]) / soargain) / npix / nmosaic[nmosaic > 0]
    noise = fits.PrimaryHDU(data = np.sqrt(varnew), header = headernew)
    noise.writeto(filename.replace(".fits", "_noise.fits"), clobber = True)

    # run sextractor on projected image
    command = "sex %s -CATALOG_NAME %s-catalogue.dat -BACK_SIZE %i -VERBOSE_TYPE QUIET" % (filename, filename, backsize)
    print command
    os.system(command)
    
    # load catalogue
    try:
        (x, y, flux, e_flux, r, flag) = np.loadtxt("%s-catalogue.dat" % (filename), usecols = (1, 2, 5, 6, 8, 9)).transpose()
    except:
        print "Sextractor failed"
        sys.exit()
    
    # load mosaic image to mask training stars
    filenmosaic = "%s/%s_%s_%02i_nmosaic_DECam_o%i.fits" % (outdir, supernova, filter, CCDSOAR, order)
    nmosaic = np.array(fits.open(filenmosaic)[0].data, dtype = int)

    # create x and y array to use nmosaic mask
    xidx = np.zeros((nx1, ny1), dtype = 'int32')
    yidx = np.zeros((ny1, nx1), dtype = 'int32')
    xidx[:] = np.arange(ny1)
    xidx = xidx.transpose()
    yidx[:] = np.arange(nx1)

    # plot reference image with stars with nmosaic < 5 masked out
    mask = (flagref <= 4)
    nref = np.empty_like(xref, dtype = int)
    pixmax1 = np.empty_like(xref)
    pixmax2 = np.empty_like(xref)
    if doplot:
        fig, ax = plt.subplots()
        ax.imshow(dataref, interpolation = 'nearest', cmap = 'gray', clim = (np.percentile(dataref, 1), np.percentile(dataref, 99)), origin = 'lower')
        ax.scatter(x[flag <= 4], y[flag <= 4], marker = 'o', c = 'b', facecolors = 'none', lw = 0.2, s = 20)
        ax.scatter(xref[mask], yref[mask], marker = 'o', c = 'r', s = 5, facecolors = 'none', lw = 0.1)

    # initialize kernel
    psf1s = None
    psf2s = None
    f1sel = []
    e_f1sel = []
    f2sel = []
    e_f2sel = []
    r1sel = []
    r2sel = []

    # update mask
    for isource in range(len(xref)):
        nref[isource] = nmosaic[int(yref[isource]) - 1, int(xref[isource]) - 1]
        pixmax1[isource] = dataref[int(yref[isource]) - 1, int(xref[isource]) - 1]
        pixmax2[isource] = datanew[int(yref[isource]) - 1, int(xref[isource]) - 1]
    mask = np.array(mask & (nref >= np.max(nmosaic) / 2.))
    maskfinal = np.zeros(np.shape(mask), dtype = bool)

    if doplot:
        ax.scatter(xref[mask], yref[mask], marker = 'o', c = 'b', s = 5, facecolors = 'none', lw = 0.5)
        plt.savefig("test_xrefyref_o%i.png" % order, dpi = 300)

    # select stars to train kernel
    for isource in range(len(xref)):

        if not mask[isource] or pixmax1[isource] > 20000 or pixmax2[isource] > 20000:
            continue

        # remove non isolated stars
        distall = np.sqrt((xref[isource] - xref[mask & (xref != xref[isource]) & (yref != yref[isource])])**2 + (yref[isource] - yref[mask & (xref != xref[isource]) & (yref != yref[isource])])**2)
        if np.min(distall) < 50:
            continue

        # remove stars close to the SN position
        print xref[isource], yref[isource], coords[1], coords[0]
        distSN = np.sqrt((xref[isource] - coords[1])**2 + (yref[isource] - coords[0])**2)
        if distSN < 50:
            print "Star too close to the SN"
            continue
        
        # take stamps and check size
        psf1 = dataref[int(yref[isource]) - dn: int(yref[isource]) + dn, int(xref[isource]) - dn: int(xref[isource]) + dn]
        psf2 = datanew[int(yref[isource]) - dn: int(yref[isource]) + dn, int(xref[isource]) - dn: int(xref[isource]) + dn]
        if np.shape(psf1) != (2 * dn, 2 * dn) or np.shape(psf2) != (2 * dn, 2 * dn):
            continue

        # remove stars near the edges of the image
        if np.median(psf1) == 0 or np.median(psf2) == 0:
            print "It appears that the image is close to the edge of one of the cameras"
            continue

        # remove bright stars on the edges of the stamp or close to the supernova position
        if rs2Dstars.flatten()[np.argmax(psf1.flatten())] >= 5 or rs2Dstars.flatten()[np.argmax(psf2.flatten())] >= 5:
            continue

        # find best match in new image
        distcat = np.sqrt((xref[isource] - x)**2 + (yref[isource] - y)**2)
        bestmatch = np.argmin(distcat)
        if distcat[bestmatch] > 4:
            continue

        # save fluxes
        f1sel.append(fluxref[isource])
        e_f1sel.append(e_fluxref[isource])
        f2sel.append(flux[bestmatch])
        e_f2sel.append(e_flux[bestmatch])
        r1sel.append(rref[isource])
        r2sel.append(r[bestmatch])
        
        # print summary
        #print isource, fluxref[isource], e_fluxref[isource], pixmax1[isource], pixmax2[isource], fluxref[isource], flux[bestmatch], rref[isource], r[bestmatch]

        # plot stars
        fig, ax = plt.subplots(2, 1)
        ax[0].imshow(psf1, interpolation = 'nearest', clim = (np.percentile(psf1, 1), np.percentile(psf1, 99)), origin = 'lower')
        ax[1].imshow(psf2, interpolation = 'nearest', clim = (np.percentile(psf2, 1), np.percentile(psf2, 99)), origin = 'lower')
        plt.savefig(filename.replace(".fits", "_kernelstar_%05i_o%i.png" % (isource, order)))

        # save psfs
        if psf1s == None:
            psf1s = psf1
            psf2s = psf1
        else:
            psf1s = np.dstack([psf1s, psf1])
            psf2s = np.dstack([psf2s, psf2])

    # count selected stars
    nstars = np.shape(psf1s)[2]
    print "Number of selected stars: %i" % nstars

    # select only stars for empirical psf
    r1sel = np.array(r1sel)
    r2sel = np.array(r2sel)
    f1sel = np.array(f1sel)
    f2sel = np.array(f2sel)
    ravg = (r1sel + r2sel) / 2.
    maskrs = (ravg < np.percentile(ravg, 75)) & (f1sel > 0) & (f2sel > 0)
        
    # decide which image to convolve
    if np.median(r1sel[maskrs]) <= np.median(r2sel[maskrs]):
        conv1st = True
    else:
        conv1st = False
    conv1st = True
    print "conv1st:", conv1st

    # infere empirical psf
    if conv1st:
        psfs = np.array(psf2s)
    else:
        psfs = np.array(psf1s)
    for i in range(nstars):
        if maskrs[i]:
            psfs[:, :, i] = psfs[:, :, i] / np.abs(np.sum(psfs[:, :, i]))
    psf = np.average(psfs[:, :, maskrs], axis = 2)
    psf[rs2Dstars > 5] = 0
    psf = psf / np.sum(psf)

    # plot psf
    if doplot:
        fig, ax = plt.subplots()
        im = ax.imshow(psf, interpolation = 'nearest', origin = 'lower')
        fig.colorbar(im)
        plt.savefig("psf_o%i.png" % order)

    # start building kernel
    X = np.zeros((nstars * npsf2, nvar))
    Y = np.zeros(nstars * npsf2)

    # loop among stars
    for i in range(nstars):

        psf1 = psf1s[:, :, i]
        psf2 = psf2s[:, :, i]

        # fill kernel equations
        for k in range(nf):
            for l in range(nf):
                ivar = ivarf[k, l]
                if ivar == -1:
                    continue
                if conv1st:
                    X[i * npsf2: (i + 1) * npsf2, ivar] \
                        = X[i * npsf2: (i + 1) * npsf2, ivar] + psf1[ieq2i + k, ieq2j + l]
                else:
                    X[i * npsf2: (i + 1) * npsf2, ivar] \
                        = X[i * npsf2: (i + 1) * npsf2, ivar] + psf2[ieq2i + k, ieq2j + l]
        if conv1st:
            Y[i * npsf2: (i + 1) * npsf2] = psf2[nfh: -(nfh + 1), nfh: -(nfh + 1)][ieq2i, ieq2j]
        else:
            Y[i * npsf2: (i + 1) * npsf2] = psf1[nfh: -(nfh + 1), nfh: -(nfh + 1)][ieq2i, ieq2j]


    # solve filter
    mat = np.dot(X.transpose(), X)
    rhs = np.dot(X.transpose(), Y)
    solvars = scipylinalg.solve(mat, rhs)

    # recover filter
    solfilter = np.zeros((nf, nf))
    for k in range(nf):
        for l in range(nf):
            ivar = int(ivarf[k, l])
            if ivar == -1:
                solfilter[k, l] = 0
                continue
            solfilter[k, l] = solvars[ivar]
    normfilter = np.sum(solfilter.flatten())
    print "normfilter: ", normfilter
    
    # plot filter
    if doplot:
        fig, ax = plt.subplots()
        ax.imshow(solfilter, interpolation = 'nearest', origin = 'lower')
        plt.savefig("solfilter_o%i.png" % order)

    # Compute normalization factor
    fig, ax = plt.subplots()
    ax.errorbar(f1sel, f2sel, xerr = e_f1sel, yerr = e_f2sel, marker = '.', lw = 0, elinewidth = 1)
    normcomp = np.sum(f1sel * f2sel) / np.sum(f1sel * f1sel)
    ax.plot([0, np.max(f1sel)], [0, np.max(f1sel) * normcomp], ls = ':')
    if conv1st:
        ax.plot([0, np.max(f1sel)], [0, np.max(f1sel) * normfilter])
    else:
        ax.plot([0, np.max(f1sel) * normfilter], [0, np.max(f1sel)])
#    solfilter = solfilter / normcomp * normfilter
    ax.set_xlabel("DECAM [ADU]")
    ax.set_ylabel("SOAR [ADU]")
    plt.savefig(filename.replace(".fits", "_fluxcalib_o%i.png" % order))

    # prepare to apply kernel
    npartx = 1
    nparty = 1
    dxconv = (nx1 - nf) / nparty
    dyconv = (ny1 - nf) / npartx

    ipart = 0
    jpart = 0
    i1 = nfh + dxconv * ipart
    i2 = i1 + dxconv
    j1 = nfh + dyconv * jpart
    j2 = j1 + dyconv

    # make sure arrays are in float32
    dataref = float32(dataref)
    datanew = float32(datanew)
    
    # apply kernel
    print "Doing convolution..."
    if conv1st:
        conv2fast.conv(dxconv + nf - 1, dyconv + nf - 1, dxconv, dyconv, nf, solfilter, dataref[i1 - nfh: i2 + nfh, j1 - nfh: j2 + nfh])
    else:
        conv2fast.conv(dxconv + nf - 1, dyconv + nf - 1, dxconv, dyconv, nf, solfilter, datanew[i1 - nfh: i2 + nfh, j1 - nfh: j2 + nfh])
        
    datat = np.empty_like(datanew)

    # recover and save solution
    print "Saving solution"
    datat[i1: i2, j1: j2] = conv2fast.iout[0: dxconv, 0: dyconv]
    conv = fits.PrimaryHDU(data = float32(datat), header = headerref)
    conv.writeto(filename.replace(".fits", "_conv.fits"), clobber = True)

    # difference
    if conv1st:
        diff = datanew - datat
    else:
        diff = datat - dataref
    difffits = fits.PrimaryHDU(data = float32(diff), header = headerref)
    difffits.writeto(filename.replace(".fits", "_diff.fits"), clobber = True)
    
    # variance of the difference
    if conv1st:
        vardiff = varref * normfilter**2 + varnew
    else:
        vardiff = varref + varnew * normfilter**2
    
    # variance
    varfits = fits.PrimaryHDU(data = float32(vardiff), header = headerref)
    varfits.writeto(filename.replace(".fits", "_vardiff.fits"), clobber = True)

    # snr
    snrfits = fits.PrimaryHDU(data = float32(diff / np.sqrt(vardiff)), header = headerref)
    snrfits.writeto(filename.replace(".fits", "_snrdiff.fits"), clobber = True)

    # save psf and conv1st
    np.save(filename.replace(".fits", "_conv1st.npy"), conv1st)
    np.save(filename.replace(".fits", "_psf.npy"), psf)
    np.save(filename.replace(".fits", "_solfilter.npy"), solfilter)
    
if dophotometry:

    # filename
    filename = "%s/%s_%s_%02i_nosky_DECam_o%i.fits" % (outdir, supernova, filter, CCDSOAR, order)
    MJD = float(fits.open(filename)[0].header['MJD-OBS'])

    # open metadata
    conv1st = np.load(filename.replace(".fits", "_conv1st.npy"))
    psf = np.load(filename.replace(".fits", "_psf.npy"))
    solfilter = np.load(filename.replace(".fits", "_solfilter.npy"))
    normpsf = np.sum(psf)
    normfilter = np.sum(solfilter)
    print conv1st, normfilter, normpsf

    # plot psf and solfilter
    if doplot:
        fig, ax = plt.subplots()
        im = ax.imshow(psf, interpolation = 'nearest', origin = 'lower')
        fig.colorbar(im)
        plt.savefig("psf_o%i.png" % order)
        
        fig, ax = plt.subplots()
        im = ax.imshow(solfilter, interpolation = 'nearest', origin = 'lower')
        fig.colorbar(im)
        plt.savefig("solfilter_o%i.png" % order)

    print rs2Dstars.flatten()[np.argmax(psf.flatten())]

    # open original image, convolved image and subtraction image
    imref = fits.open(fitsref)[0].data
    imnew = fits.open(filename)[0].data
    imconv = fits.open(filename.replace(".fits", "_conv.fits"))[0].data

    # open difference and variance maps
    diff = fits.open(filename.replace(".fits", "_diff.fits"))[0].data
    var = fits.open(filename.replace(".fits", "_vardiff.fits"))[0].data
    snr = diff / var

    # do better centering (maybe try doing some previous filtering?)
    (ix, iy) = coords
    aux = ndimage.filters.gaussian_filter(diff[ix - 2: ix + 2, iy - 2: iy + 2], 4.)
    (dx, dy) = np.unravel_index(np.argmax(aux), np.shape(aux))
    print "dx: %s, dy: %s" % (dx, dy)
    ix = ix + (dx - 1) * 0.5
    iy = iy + (dy - 1) * 0.5
    
    # mask for region to be included in the photometry
    maskphoto = (psf >= 0.2 * np.max(psf)) & (rs2Dstars < 4)
    imref = imref[ix - dn: ix + dn, iy - dn: iy + dn]
    imnew = imnew[ix - dn: ix + dn, iy - dn: iy + dn]
    imconv = imconv[ix - dn: ix + dn, iy - dn: iy + dn]
    diff = diff[ix - dn: ix + dn, iy - dn: iy + dn]
    var = var[ix - dn: ix + dn, iy - dn: iy + dn]
    snr = snr[ix - dn: ix + dn, iy - dn: iy + dn]

    # raw sum
    print "Raw sum:", np.sum(diff[maskphoto])

    # compute weights
    weights = np.array(psf / var / np.sum(psf[maskphoto]**2 / var[maskphoto]))
    weights[np.invert(maskphoto)] = 0
    
    # psf masked
    psfmasked = np.array(psf)
    psfmasked[np.invert(maskphoto)] = 0

    # plot regions around SNe
    fig, ax = plt.subplots(2, 4, figsize = (12, 6))
    im = ax[0, 0].imshow(imref, interpolation = 'nearest', origin = 'lower')
    ax[0, 0].set_title("imref", fontsize = 6)
    im = ax[0, 1].imshow(imnew, interpolation = 'nearest', origin = 'lower')
    ax[0, 1].set_title("imnew", fontsize = 6)
    if conv1st:
        convlabel = 'imref_conv'
    else:
        convlabel = 'imnew_conv'
    im = ax[0, 2].imshow(imconv, interpolation = 'nearest', origin = 'lower')
    ax[0, 2].set_title(convlabel, fontsize = 6)
    im = ax[0, 3].imshow(diff, interpolation = 'nearest', origin = 'lower')
    ax[0, 3].contour(psfmasked, alpha = 0.3)
    ax[0, 3].set_title("diff", fontsize = 6)
    im = ax[1, 0].imshow(weights, interpolation = 'nearest', origin = 'lower')
    ax[1, 0].contour(psfmasked)
    ax[1, 0].set_title("weights", fontsize = 6)
    im = ax[1, 1].imshow(psfmasked, interpolation = 'nearest', origin = 'lower')
    ax[1, 1].contour(psfmasked)
    ax[1, 1].set_title("psf", fontsize = 6)
    im = ax[1, 2].imshow(var, interpolation = 'nearest', origin = 'lower')
    ax[1, 2].contour(rs2Dstars)
    ax[1, 2].set_title("variance", fontsize = 6)
    im = ax[1, 3].imshow(snr, interpolation = 'nearest', origin = 'lower')
    ax[1, 3].contour(psfmasked)
    ax[1, 3].set_title("snr", fontsize = 6)
    plt.savefig(filename.replace(".fits", "_supernova.png"), bbox_inches = 'tight', pad_inches = 0.01)

#    # plot weights
#    fig, ax = plt.subplots()
#    im = ax.imshow(weights, interpolation = 'nearest', origin = 'lower')
#    fig.colorbar(im)
#    plt.savefig("weights.png")
#
#    # plot variance 
#    var[np.invert(maskphoto)] = 0
#    fig, ax = plt.subplots()
#    im = ax.imshow(var, interpolation = 'nearest', origin = 'lower', clim = (0, 300))
#    fig.colorbar(im)
#    plt.savefig("var.png")
#
#    # plot psf
#    psf[np.invert(maskphoto)] = 0
#    fig, ax = plt.subplots()
#    im = ax.imshow(psf, interpolation = 'nearest', origin = 'lower')
#    fig.colorbar(im)
#    plt.savefig("psf.png")

    # compute flux
    flux = np.sum(weights * diff)
    e_flux = np.sqrt(np.sum(weights**2 * var))
    print "Weighted flux: %s +- %s (SNR: %s)" % (flux, e_flux, flux / e_flux)
    if conv1st:
        flux = flux / normfilter
        e_flux = e_flux / normfilter
    print "Weighted flux, normed: %s +- %s (SNR: %s)" % (flux, e_flux, flux / e_flux)

    # compute flux
    if doplot:
        fig, ax = plt.subplots()
        im = ax.imshow(weights * diff, interpolation = 'nearest', origin = 'lower')
        fig.colorbar(im)
        plt.savefig("flux-weights_o%i.png" % order)
        diff[np.invert(maskphoto)] = 0
        fig, ax = plt.subplots()
        im = ax.imshow(diff, interpolation = 'nearest', origin = 'lower', clim = (0, 200))
        fig.colorbar(im)
        plt.savefig("flux_o%i.png" % order)

    # open CCD numbers file
    CCDn = {}
    (CCDstring, CCDnumber) = np.loadtxt("CCDnumbers.dat", dtype = str).transpose()
    CCDnumber = np.array(CCDnumber, dtype = int)
    for i in range(len(CCDstring)):
        CCDn[CCDstring[i]] = CCDnumber[i]

    # open zero point file
    (ID, filter, a, aerr, b, berr, k, kerr) = np.loadtxt("zeropoints_g.txt", dtype = str).transpose()
    ID = np.array(ID, dtype = int)
    a = np.array(a, dtype = float)
    aerr = np.array(aerr, dtype = float)
    b = np.array(b, dtype = float)
    berr = np.array(berr, dtype = float)
    k = np.array(k, dtype = float)
    kerr = np.array(kerr, dtype = float)

    # exposure time and airmass
    exptimeref = float(headerref['EXPTIME'])
    airmassref = float(headerref['AIRMASS'])
    print exptimeref, airmassref

    # convert flux to magnitudes
    idx = CCDn[CCD] - 1
    mag = np.array(-2.5 * np.log10(flux) + 2.5 * np.log10(exptimeref) - a[idx] - k[idx] * airmassref)
    if flux > e_flux:
        mag_1 = np.array(-2.5 * np.log10(flux - e_flux) + 2.5 * np.log10(exptimeref) - a[idx] - k[idx] * airmassref)
        e_mag_1 = mag_1 - mag
    else:
        mag_1 = 0
        e_mag_1 = 0
    mag_2 = np.array(-2.5 * np.log10(flux + e_flux) + 2.5 * np.log10(exptimeref) - a[idx] - k[idx] * airmassref)
    e_mag_2 = mag - mag_2
    print "Weighted flux (mag): %s (-%s +%s)" % (mag, e_mag_2, e_mag_1)

    # save results
    np.savetxt(filename.replace(".fits", "_flux.dat"), [MJDref, MJD, MJD - MJDref, flux, e_flux, mag, e_mag_2, e_mag_1])
