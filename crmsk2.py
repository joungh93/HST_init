#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 8 13:51:43 2020

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
from astropy.nddata import Cutout2D
from astropy import wcs
from astropy.stats import sigma_clipped_stats
from tqdm import trange


# ----- Make an image list ----- #
alw = glob.glob('alw*_drc.fits')
alw = sorted(alw)
diT = "Trim/"
objname = "JFG2"
jfg = glob.glob(diT+'i*.fits')
jfg = sorted(jfg)


# ----- Circular masking ----- #
msk = [(186., 124., 4.0),
       (190., 119., 4.0),
       (207., 127., 3.5),
       (218., 127., 4.0),
       (58., 69., 3.0),
       (141., 55., 3.0),
       (187., 56., 5.0),
       (54., 22., 3.0)]
       # (x, y, r) in jfg
umsk = (138, 107)  # in jfg

kernel = Gaussian2DKernel(5)

for i in np.arange(len(jfg)):
	filt = jfg[i].split('/')[1].split('.fits')[0].split('_')[1]

	print("Masking "+jfg[i]+"...")

	img, hdr = fits.getdata(jfg[i], header=True, ext=0)
	w = wcs.WCS(hdr)

	img2 = copy.deepcopy(img)
	fl_msk = np.zeros((img.shape[0], img.shape[1]))

	x = np.arange(0, img.shape[1], 1)
	y = np.arange(0, img.shape[0], 1)
	xx, yy = np.meshgrid(x, y, sparse=True)

	# Circular masking
	for j in trange(len(msk)):
		z = (xx-(msk[j][0]-1))**2.0 + (yy-(msk[j][1]-1))**2.0 - (1.0*msk[j][2])**2.0
		mskpix = (z <= 0.0)
		img2[mskpix] = np.nan
		fl_msk[mskpix] = 1

	img_conv = convolve(img2, kernel)
	img2[fl_msk == 1] = img_conv[fl_msk == 1]

	# Unmasking
	img_alw, hd1_alw = fits.getdata(alw[i], header=True, ext=1)
	w_alw = wcs.WCS(hd1_alw)
	ra_umsk, dec_umsk = w.wcs_pix2world(umsk[0], umsk[1], 1)
	x_alw, y_alw = w_alw.wcs_world2pix(ra_umsk, dec_umsk, 1)

	ch = Cutout2D(data=img_alw, position=(x_alw, y_alw), size=(21, 31))
	img2[umsk[1]-1-int(21/2):umsk[1]-1+int(21/2)+1,
	     umsk[0]-1-int(31/2):umsk[0]-1+int(31/2)+1] = ch.data

	# Region masking
	avg, med, std = sigma_clipped_stats(img2, sigma=3.0, maxiters=10)
	np.random.seed(0)
	img2[np.isnan(img2) == True] = np.random.normal(med, std, img2[np.isnan(img2) == True].shape)
	img2[:20, :] = np.random.normal(med, std, img2[:20, :].shape)
	img2[:60, 181:] = np.random.normal(med, std, img2[:60, 181:].shape)
	img2[50:100, 250:300] = np.random.normal(med, std, img2[50:100, 250:300].shape)

	# Creating new images
	fits.writeto(diT+'mi'+objname+'_'+filt+'.fits', img2, hdr, overwrite=True) 


# Printing the running time  
print("--- %s seconds ---" % (time.time() - start_time))