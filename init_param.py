#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 14:30:58 2020

@author: jlee
"""


import numpy as numpy
import glob, os


# ----- Initial science images ----- #


# mk_fits.py
sci_dir = 'Images/'
sci_img = glob.glob(sci_dir+'*.fits')
sci_img = sorted(sci_img)
flt_Ref = '606'
num_Ref = '1'

# apply_sep0.py
tmp_dir = 'temp/'

# mk_comb.py
pixscl = 0.05
mch_tol = 2.0
min_npoi = 10
