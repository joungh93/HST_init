#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 8 17:59:16 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
import pandas as pd
import copy
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import wcs


# ----- Initial settings ----- #
diT = "Trim/"
objname = "JFG2"


# ----- Make an image list ----- #
m1alw = glob.glob('m1alw*_drc.fits')
m1alw = sorted(m1alw)


# ----- Making cutout images ----- #
for i in np.arange(len(m1alw)):
	filt = m1alw[i].split('_drc.fits')[0].split('m1alw')[1]

	# hd0 = fits.getheader(m1alw[i], ext=0)
	img, hd1 = fits.getdata(m1alw[i], header=True, ext=1)
	w = wcs.WCS(hd1)

	chdu = Cutout2D(data=img, position=(430, 150), size=(300, 300), wcs=w)
	fits.writeto(diT+'i'+objname+'_'+filt+'.fits', chdu.data, chdu.wcs.to_header(),
		         overwrite=True)

	# nhd0 = fits.PrimaryHDU()
	# nhd0.header = hd0

	# nhd1 = fits.ImageHDU()
	# nhd1.data = chdu.data
	# nhd1.header = chdu.wcs.to_header()
	# for keywords in ['XTENSION', 'PCOUNT', 'GCOUNT']:
	# 	nhd1.header[keywords] = hd1[keywords]

	# nhdu = fits.HDUList([nhd0, nhd1])
	# nhdu.writeto(diT+objname+'_'+filt+'.fits', overwrite=True) 


# Printing the running time  
print("--- %s seconds ---" % (time.time() - start_time))
