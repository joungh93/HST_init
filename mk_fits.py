#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:41:35 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
from astropy.io import fits
import glob, os
from astroscrappy import detect_cosmics
from astropy.stats import sigma_clipped_stats
from reproject import reproject_interp

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

import init_param as ip


# ----- Making the initial science images ----- #
os.system('rm -rfv ./*.fits ./*.reg ./*.csv')    # Initialization
os.system('rm -rfv '+ip.tmp_dir)
os.system('mkdir '+ip.tmp_dir)

filt = np.array([])
for i in np.arange(len(ip.sci_img)):
    ### Reading data
    hdu = fits.open(ip.sci_img[i])
    dat_sci= hdu[0].data
    hdr_sci = hdu[0].header
    if (hdr_sci['INSTRUME'] == 'ACS'):
        if (hdr_sci['FILTER1'][0:5] != 'CLEAR'):
            hfilt = hdr_sci['FILTER1']
        if (hdr_sci['FILTER2'][0:5] != 'CLEAR'):
            hfilt = hdr_sci['FILTER2']
    if (hdr_sci['INSTRUME'] == 'WFC3'):
        hfilt = hdr_sci['FILTER']
    filt = np.append(filt, hfilt)
    img_already_written = glob.glob('./'+hfilt[1:4]+'_*.fits')
    itag_num = len(img_already_written)+1
    new_img_name = hfilt[1:4]+'_'+str(itag_num)+'.fits'

    ### Cosmic ray rejection

    # GAIN & Readout noise
    gain_value = []
    for k in ['ATODGNA', 'ATODGNB', 'ATODGNC', 'ATODGND']:
        gain_value.append(hdr_sci[k])
    epadu = np.mean(np.array(gain_value))
    
    read_value = []
    for k in ['READNSEA', 'READNSEB', 'READNSEC', 'READNSED']:
        read_value.append(hdr_sci[k])
    rdnoi = np.mean(np.array(read_value))
    exptime = hdr_sci['EXPTIME']

    # Creating a masking array
    msk = np.zeros_like(dat_sci, dtype=bool)
    msk[np.isnan(dat_sci) == True] = True

    avg, med, std = sigma_clipped_stats(dat_sci, mask=msk, sigma=3.0,
                                        maxiters=5, std_ddof=1)
    if (med-avg < 0.3*std):
        msky = med
    else:
        msky = 2.5*med - 1.5*avg
    msk[dat_sci < msky-1.0*std] = True

    # Cleaning the images
    print('Cleaning '+new_img_name+' ...')
    LACOSMIC_KEYS = dict(sigclip=4.5, sigfrac=0.3, objlim=5.0,
                         satlevel=65535.0, pssl=0.0, niter=4,
                         cleantype='medmask', fsmode='median', psfmodel='gauss',
                         psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765)
    crmask, cleanarr = detect_cosmics(exptime*dat_sci, inmask=msk,
                                      gain=epadu, readnoise=rdnoi, **LACOSMIC_KEYS)
    cleanarr /= epadu*exptime
    print(f'{np.sum(crmask):d} cosmic ray pixels removed\n')

    ### Writing the output image
    fits.writeto(new_img_name, cleanarr, hdr_sci, overwrite=True)

filt = np.unique(filt)
# Output: [FILTER]_[NUMBER].fits


# ----- Registering WCS for each image ----- #
hdu_ref = fits.open(ip.flt_Ref+'_'+ip.num_Ref+'.fits')[0]

for i in np.arange(len(filt)):
    flt_num = filt[i][1:4]
    img = glob.glob(flt_num+'_*.fits')
    img = sorted(img)
    for j in np.arange(len(img)):
        print('Reprojecting & writing '+img[j]+' ...')
        hdu = fits.open(img[j])[0]
        if (img != ip.flt_Ref+'_'+ip.num_Ref+'.fits'):
            array, footprint = reproject_interp(hdu, hdu_ref.header)
            for k in ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                      'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
                hdu.header[k] = hdu_ref.header[k]
            fits.writeto('w'+flt_num+'_{0:d}.fits'.format(j+1),
                         array, hdu.header, overwrite=True)
        else:
            os.system('cp -rpv '+ip.flt_Ref+'_'+ip.num_Ref+'.fits w'+ip.flt_Ref+'_'+ip.num_Ref+'.fits')
# Output: w[FILTER]_[NUMBER].fits

# for i in np.arange(len(filt)):
#     flt_num = filt[i][1:4]
#     os.system('rm -rfv '+flt_num+'_*.fits')


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
