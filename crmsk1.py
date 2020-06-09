#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 17:16:10 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
import pandas as pd
import copy
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from tqdm import trange


# ----- Make an image list ----- #
alw = glob.glob('alw*_drc.fits')
alw = sorted(alw)


# ----- CR pixels ----- #
dir_CR = '/data/jlee/DATA/HLA/McPartland+16/MACS1752/Img_NE/Regions/'
cr = np.genfromtxt(dir_CR+'cos2.reg', dtype=None, encoding='ascii', names=('x','y','num'))


# ----- Loading the pickle files ----- #
sep1_pkl = glob.glob('d_sep1_*.pkl')
sep1_pkl = sorted(sep1_pkl)
for i in np.arange(len(sep1_pkl)):
    dat = pd.read_pickle(sep1_pkl[i])
    exec(sep1_pkl[i].split('.')[0]+' = dat')


# ----- Masking images ----- #
kernel = Gaussian2DKernel(5)

for i in np.arange(len(alw)):

    print("Masking "+alw[i]+"...")

    exec('ds = '+sep1_pkl[i].split('.')[0])

    hd0 = fits.getheader(alw[i], ext=0)
    img, hd1 = fits.getdata(alw[i], header=True, ext=1)

    img2 = copy.deepcopy(img)
    fl_msk = np.zeros((img.shape[0], img.shape[1]))

    x = np.arange(0, img.shape[1], 1)
    y = np.arange(0, img.shape[0], 1)
    xx, yy = np.meshgrid(x, y, sparse=True)

    for j in trange(len(cr)):
        idx = cr['num'][j]-1
        z = ds.iloc[idx]['cxx']*(xx-(cr['x'][j]-1))**2.0 + \
            ds.iloc[idx]['cyy']*(yy-(cr['y'][j]-1))**2.0 + \
            ds.iloc[idx]['cxy']*(xx-(cr['x'][j]-1))*(yy-(cr['y'][j]-1))- \
            (1.0*ds.iloc[idx]['kron'])**2.0
        mskpix = (z <= 0.0)
        img2[mskpix] = np.nan
        fl_msk[mskpix] = 1
    img_conv = convolve(img2, kernel)
    img2[fl_msk == 1] = img_conv[fl_msk == 1]

    nhd0 = fits.PrimaryHDU()
    nhd0.header = hd0

    nhd1 = fits.ImageHDU()
    nhd1.data = img2
    nhd1.header = hd1

    nhdu = fits.HDUList([nhd0, nhd1])
    nhdu.writeto('m1'+alw[i], overwrite=True)   


# Printing the running time  
print("--- %s seconds ---" % (time.time() - start_time))