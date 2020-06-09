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
from reproject import reproject_interp
import astroalign as aa


# ----- Raw images ----- #
raw_dir = 'Raw/'
raw_image = glob.glob("Raw/*.fits")
raw_image = sorted(raw_image)


# ----- Changing image name w/ filter information from the image headers ----- #
filt = []
for i in np.arange(len(raw_image)):
    hdu = fits.open(raw_image[i])
    hdr = hdu[0].header

    if (hdr['INSTRUME'] == 'ACS'):
        if (hdr['FILTER1'][0:5] != 'CLEAR'):
            hfilt = hdr['FILTER1']
        if (hdr['FILTER2'][0:5] != 'CLEAR'):
            hfilt = hdr['FILTER2']
    if (hdr['INSTRUME'] == 'WFC3'):
        hfilt = hdr['FILTER']

    filt.append(hfilt)

    dat_sci = hdu[1].data
    hdr_sci = hdu[1].header
    hdr_sci.remove('EXTNAME')
    hdr_sci.remove('EXTVER')

    nhd0 = fits.PrimaryHDU()
    nhd0.header = hdr

    nhd1 = fits.ImageHDU()
    nhd1.data = dat_sci
    nhd1.header = hdr_sci

    nhdu = fits.HDUList([nhd0, nhd1])
    nhdu.writeto(hfilt[1:-1]+'_drc.fits', overwrite=True)
# Output: *_drc.fits


# ----- Registering WCS for each image ----- #
flt_Ref = '606'
for i in np.arange(len(filt)):
	flt_num = filt[i][1:-1]
	exec("hd"+flt_num+" = fits.open('"+flt_num+"_drc.fits')[1]")
	if (flt_num == flt_Ref):
		exec("head_ref = hd"+flt_num+".header")

for i in np.arange(len(filt)):
	flt_num = filt[i][1:-1]
	hdr = fits.getheader(flt_num+'_drc.fits', ext=0)

	if (flt_num != flt_Ref):
		exec("hd_input = hd"+flt_num)
		array, footprint = reproject_interp(hd_input, head_ref)

		for j in ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2']:
			hd_input.header[j] = head_ref[j]

		nhd0 = fits.PrimaryHDU()
		nhd0.header = hdr

		nhd1 = fits.ImageHDU()
		nhd1.data = array
		nhd1.header = hd_input.header

		nhdu = fits.HDUList([nhd0, nhd1])
		nhdu.writeto('w'+flt_num+'_drc.fits', overwrite=True)

	else:
		os.system('cp -rpv '+flt_num+'_drc.fits w'+flt_num+'_drc.fits')
# Output: w*_drc.fits


# ----- Spatial alignment ----- #
flt_Ref = '606'
for i in np.arange(len(filt)):
	flt_num = filt[i][1:-1]
	exec("hd"+flt_num+" = fits.open('w"+flt_num+"_drc.fits')[1]")
	if (flt_num == flt_Ref):
		exec("hd_ref = hd"+flt_num)

for i in np.arange(len(filt)):
	flt_num = filt[i][1:-1]
	hdr = fits.getheader('w'+flt_num+'_drc.fits', ext=0)
	if (flt_num != flt_Ref):
		exec("hd_input = hd"+flt_num)
		transf, (src_lis, tar_lis) = aa.find_transform(hd_input.data, hd_ref.data)
		exec("transf_"+flt_num+" = transf")
		exec("src_"+flt_num+" = src_lis")
		exec("tar_"+flt_num+" = tar_lis")

		# Printing out the results
		exec("p = transf_"+flt_num)
		print("\n# ----- Parameters of the transformation: "+filt[i]+" band ----- #")
		print("Rotation: {:.2f} degrees".format(p.rotation * 180.0/np.pi))
		print("Scale factor: {:.2f}".format(p.scale))
		print("Translation: (x, y) = ({:.2f}, {:.2f})".format(*p.translation))
		print("Tranformation matrix:\n{}".format(p.params))

		# Writing region files 
		f = open('source_'+filt[i]+'.reg', 'w')
		for j in np.arange(src_lis.shape[0]):
			f.write('{0:.3f}  {1:.3f}'.format(src_lis[j][0], src_lis[j][1])+'\n')
		f.close()	

		# Alignment & saving new FITS files
		img_aligned, footprint = aa.register(hd_input.data, hd_ref.data)

		nhd0 = fits.PrimaryHDU()
		nhd0.header = hdr

		nhd1 = fits.ImageHDU()
		nhd1.data = img_aligned
		nhd1.header = hd_input.header

		nhdu = fits.HDUList([nhd0, nhd1])
		nhdu.writeto('alw'+flt_num+'_drc.fits', overwrite=True)

	else:
		os.system('cp -rpv w'+flt_num+'_drc.fits alw'+flt_num+'_drc.fits')
		print("\n# ----- The reference filter: "+filt[i]+" band ----- #")



# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
