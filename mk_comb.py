#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 11:57:21 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
from astropy.io import fits
import glob, os


# ----- Making a detection image ----- #
alw = glob.glob('alw*_drc.fits')
alw = sorted(alw)

dat = fits.getdata(alw[0], ext=0)
cdat = np.zeros((dat.shape[0], dat.shape[1], len(alw)))

out = 'calw_dtc.fits'
for i in np.arange(len(alw)):
	dat, hdr = fits.getdata(alw[i], header=True, ext=0)
	cdat[:,:,i] = dat
cdat2 = np.mean(cdat, axis=2)
fits.writeto(out, cdat2, hdr, overwrite=True)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
