#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 14:30:58 2020

@author: jlee
"""


import numpy as numpy
import glob, os


# ----- Initial science images ----- #

# Data structure & initial science images
sci_dir = 'Images/'    # The directory including the drizzled images
sci_img = glob.glob(sci_dir+'*.fits')    # The drizzled images
sci_img = sorted(sci_img)
tmp_dir = 'temp/'    # The directory for saving source catalogs and regions

# The reference image
flt_Ref = '606'    # The reference filter for reprojecting & aligning
num_Ref = '1'    # The image number of the reference image (Maybe you should check it visually!)

# Source matching parameters for aligning images
pixscl = 0.05    # Pixel scale: arcsec/pix
mch_tol = 2.0    # Matching tolerance [pix]
min_npoi = 10    # Minimum point sources for applying transformation
