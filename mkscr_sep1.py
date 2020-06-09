#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 13:45:24 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
from astropy.io import fits
import pandas as pd


# ----- Parameter setting ----- #
DETECT_THRESH = 2.0
DETECT_MINAREA = 3
DEBLEND_NTHRESH = 32
DEBLEND_MINCONT = 0.005
PHOT_APERTURES = [3.0, 6.0, 8.0]
# CHECKIMAGE_TYPE = ['BACKGROUND', 'FILTERED', 'APERTURES', 'SEGMENTATION']
# CHECKIMAGE_NAME = ['bgr', 'flt', 'apr', 'seg']
# ----------------------------- #


# ----- Aligned images (w/ WCS information) ----- #
alw = glob.glob('alw*_drc.fits')
alw = sorted(alw)
dtc = 'calw_dtc.fits'


# ----- Writing the script ----- #
f = open('sep1.sh','w')
f.write(' \n')
f.write('##################################'+'\n')
f.write('##### Scripts for SExtractor #####'+'\n')
f.write('##################################'+'\n')        
f.write(' \n')
for i in np.arange(len(alw)):
	flt = alw[i].split('_drc.fits')[0].split('alw')[1]
	
	hdr0 = fits.getheader(alw[i], ext=0)
	hdr1 = fits.getheader(alw[i], ext=1)

	gain = 2.0*hdr0['EXPTIME']
	zp = -2.5*np.log10(hdr1['PHOTFLAM'])-5.0*np.log10(hdr1['PHOTPLAM'])-2.408

	f.write('sex '+dtc+','+alw[i]+'[1] -c config.sex -CATALOG_NAME sep1_'+flt+'.cat ')
	f.write('-DETECT_THRESH {0:.1f} -ANALYSIS_THRESH {0:.1f} '.format(DETECT_THRESH))
	f.write('-DETECT_MINAREA {0:.1f} -DEBLEND_NTHRESH {1:d} -DEBLEND_MINCONT {2:.3f} '.format(DETECT_MINAREA, DEBLEND_NTHRESH, DEBLEND_MINCONT))
	f.write('-MAG_ZEROPOINT {0:.3f} -GAIN {1:.1f} '.format(zp, gain))
	
	photap = ''
	for j in np.arange(len(PHOT_APERTURES)):
		if (j != np.arange(len(PHOT_APERTURES))[-1]):
			photap += '{0:.1f}'.format(PHOT_APERTURES[j])+','
		else:
			photap += '{0:.1f}'.format(PHOT_APERTURES[j])
	f.write('-PHOT_APERTURES '+photap+' ')

	# check_type, check_name = '', ''
	# for j in np.arange(len(CHECKIMAGE_TYPE)):
	# 	if (j != np.arange(len(CHECKIMAGE_TYPE))[-1]):
	# 		check_type += CHECKIMAGE_TYPE[j]+','
	# 		check_name += flt+'_'+CHECKIMAGE_NAME[j]+'.fits,'
	# 	else:
	# 		check_type += CHECKIMAGE_TYPE[j]
	# 		check_name += flt+'_'+CHECKIMAGE_NAME[j]+'.fits'
	# f.write('-CHECKIMAGE_TYPE '+check_type+' -CHECKIMAGE_NAME '+check_name+' ')
	
	f.write('\n')
f.close()


# ----- Running scripts for SExtractor ----- #
os.system('sh sep1.sh')


# ----- Reading & saving the results ----- #
colname = ['x','y','num','mag_auto','merr_auto',
           'mag1','mag2','mag3','merr1','merr2','merr3',
           'kron','petro','bgr','ra','dec','cxx','cyy','cxy',
           'a','b','theta','mu0','mut','flag','fwhm','flxrad','cl']

for i in np.arange(len(alw)):
	flt = alw[i].split('_drc.fits')[0].split('alw')[1]
	dat = np.genfromtxt('sep1_'+flt+'.cat', dtype=None, comments='#', encoding='ascii',
                        names=colname)
	dat = pd.DataFrame(dat)
	dat.to_pickle('d_sep1_'+flt+'.pkl')


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))