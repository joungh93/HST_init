#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 09 16:23:57 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import astroalign as aa
import init_param as ip
from apply_sep0 import sep0


# ----- Input images: w[FITLER]_[NUMBER].fits ----- #
wimg = glob.glob('w*.fits')
wimg = sorted(wimg)


# ----- Aligning images  ----- #
s0 = sep0('w'+ip.flt_Ref+'_'+ip.num_Ref+'.fits')
df0 = s0.AutoPhot()
df0_poi = s0.Pick_PointSrc(mag_high=17.0, mag_low=27.0,
                           size_low=1.0, size_high=4.0,
                           ar_cut=0.75)
src0 = SkyCoord(ra=df0_poi['ra'].values*u.deg, dec=df0_poi['dec'].values*u.deg)

tol = ip.pixscl * ip.mch_tol / 3600.0
trans_tol = ip.mch_tol

filt = []
for i in np.arange(len(wimg)):
    wname = wimg[i].split('.fits')[0]
    if (wimg[i] != 'w'+ip.flt_Ref+'_'+ip.num_Ref+'.fits'):
        print("\nAligning "+wimg[i]+"...")

        s1 = sep0(wimg[i])
        df1 = s1.AutoPhot()
        df1_poi = s1.Pick_PointSrc(mag_high=17.0, mag_low=27.0,
                                   size_low=1.0, size_high=4.0,
                                   ar_cut=0.75)

        src1 = SkyCoord(ra=df1_poi['ra'].values*u.deg, dec=df1_poi['dec'].values*u.deg)

        idx, d2d, d3d = src1.match_to_catalog_sky(src0)
        matched = d2d.value < tol
        n_mch = np.sum(matched)
        midx1 = np.where(matched)[0]
        midx0 = idx[matched]
        print(n_mch)

        g = open(ip.tmp_dir+'poi0_mat_'+wname+'.reg','w')
        g.write('global color=red dashlist=8 3 width=2 font="helvetica 10 normal roman" ')
        g.write('select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
        g.write('fk5 \n')
        for k in np.arange(n_mch):
          g.write('circle(%.5f, %.5f, 0.5")\n' \
                  %(src0.ra.value[midx0][k], src0.dec.value[midx0][k]))
        g.close()

        xc0, yc0 = df0_poi.iloc[midx0]['x'].values, df0_poi.iloc[midx0]['y'].values
        xc1, yc1 = df1_poi.iloc[midx1]['x'].values, df1_poi.iloc[midx1]['y'].values 
        src = np.array([xc1, yc1]).T ; dst = np.array([xc0, yc0]).T
        tform = aa.estimate_transform('affine', src, dst)
        print("# ----- Parameters of the transformation ----- #")
        print(f"Rotation: {tform.rotation*180.0/np.pi:.2f} degrees")
        print("Scale factor: ({:.2f}, {:.2f})".format(*tform.scale))
        print("Translation: (x, y) = ({:.2f}, {:.2f})".format(*tform.translation))
        print(f"Tranformation matrix:\n{tform.params}")
        print("\n")

        if ((n_mch < ip.min_npoi) | (np.sqrt(np.sum(tform.translation**2.0)) > trans_tol)):
            print("Warning: please check this image visually: al"+wimg[i])
            os.system('cp -rpv '+wimg[i]+' al'+wimg[i])
        else:
            aligned_image, footprint = aa.apply_transform(tform, s1.dat, s0.dat)
            fits.writeto('al'+wimg[i], aligned_image, s1.hdr, overwrite=True)
    else:
        os.system('cp -rpv '+wimg[i]+' al'+wimg[i])
        ysize, xsize = fits.getdata(wimg[i], ext=0).shape
    filt.append(wname.split('w')[1].split('_')[0])
filt = np.array(filt)
# Output: alw[FILTER]_[NUMBER].fits


# ----- Combining images  ----- #
ufilt, nfilt = np.unique(filt, return_counts=True)

for i in np.arange(len(ufilt)):
	print("Combining "+ufilt[i]+" band...")
	cdat = np.zeros((ysize, xsize))
	cexp = 0.
	for j in np.arange(nfilt[i]):
		dat, hdr = fits.getdata('alw'+ufilt[i]+f'_{j+1:d}'+'.fits', ext=0, header=True)
		dat[np.isnan(dat) == True] = 0.0
		cdat += dat
		cexp += hdr['EXPTIME']
	cdat /= nfilt[i]
	hdr['EXPTIME'] = cexp
	fits.writeto('calw'+ufilt[i]+'.fits', cdat, hdr, overwrite=True)
# Output: calw[FILTER].fits


# ----- Creating the detection image  ----- #
calw_img = glob.glob('calw*.fits')
calw_img = sorted(calw_img)

cdat = np.zeros((ysize, xsize))
out = 'calw_dtc.fits'
for i in np.arange(len(calw_img)):
    print("Creating the detection image...")
    d, h = fits.getdata(calw_img[i], header=True, ext=0)
    if (h['INSTRUME'] == 'ACS'):
        cdat += d
fits.writeto(out, cdat, h, overwrite=True)

# os.system('rm -rfv alw*.fits')


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
