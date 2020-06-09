#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 17:16:13 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker


# ----- Loading the pickle files ----- #
sep1_pkl = glob.glob('d_sep1_*.pkl')
sep1_pkl = sorted(sep1_pkl)
for i in np.arange(len(sep1_pkl)):
	dat = pd.read_pickle(sep1_pkl[i])
	exec(sep1_pkl[i].split('.')[0]+' = dat')


# ----- Figure & grid setting ----- #
band = ['F435W','F606W','F814W']
color = ['F435W-F606W', 'F606W-F814W', 'F606W-F814W']

mag_lim, mag_tick = [17, 29], [18, 20, 22, 24, 26, 28]
merr_lim, merr_tick = [0., 0.5], [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
cidx_lim, cidx_tick = [-0.9,3.5], [-1., 0., 1., 2., 3.]
flxr_lim, flxr_tick = [-2.0,10.0], [-2., 0., 2., 4., 6., 8., 10.]
fwhm_lim, fwhm_tick = [-4.0,20.0], [-4., 0., 4., 8., 12., 16., 20.]
mu0_lim, mu0_tick = [15, 26], [16, 18, 20, 22, 24, 26]
col_lim, col_tick = [-1.0, 2.5], [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5]

for i in np.arange(len(band)):
	exec("ds = d_sep1_"+band[i][1:4])

	fig = plt.figure(i+1, figsize=(11,7))
	gs = GridSpec(2, 3, left=0.08, bottom=0.12, right=0.97, top=0.97,
	              height_ratios=[1.,1.], width_ratios=[1.,1.,1.], hspace=0.10, wspace=0.30)

	ax1 = fig.add_subplot(gs[0,0])
	ax2 = fig.add_subplot(gs[0,1])
	ax3 = fig.add_subplot(gs[1,0])
	ax4 = fig.add_subplot(gs[1,1])
	ax5 = fig.add_subplot(gs[0,2])
	ax6 = fig.add_subplot(gs[1,2])

	# Axis 1 : mag - magerr
	ax = ax1
	ax.set_xticks(mag_tick)
	ax.set_xticklabels([])
	ax.set_yticks(merr_tick)
	ax.set_yticklabels(merr_tick, fontsize=12.0)
	ax.set_ylabel(band[i]+' error', fontsize=12.0)
	ax.set_xlim(mag_lim)
	ax.set_ylim(merr_lim)
	ax.tick_params(width=1.5, length=8.0)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
	ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
	ax.tick_params(width=1.5,length=5.0,which='minor')
	for axis in ['top','bottom','left','right']:
	    ax.spines[axis].set_linewidth(1.5)

	X = ds['mag_auto'].values
	Y = ds['merr_auto'].values
	ax.plot(X, Y, 'o', ms=2.0, color='darkgray', alpha=0.7)

	mag_cnd = ((X > 30.0) | (Y > 1.0))    # for cosmic ray candidates

	# Axis 2 : mag - C index
	ax = ax2
	ax.set_xticks(mag_tick)
	ax.set_xticklabels([])
	ax.set_yticks(cidx_tick)
	ax.set_yticklabels(cidx_tick, fontsize=12.0)
	ax.set_ylabel('Concentration index', fontsize=12.0)
	ax.set_xlim(mag_lim)
	ax.set_ylim(cidx_lim)
	ax.tick_params(width=1.5, length=8.0)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
	ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
	ax.tick_params(width=1.5,length=5.0,which='minor')
	for axis in ['top','bottom','left','right']:
	    ax.spines[axis].set_linewidth(1.5)

	X = ds['mag_auto'].values
	Y = ds['mag1'].values-ds['mag3'].values
	ax.plot(X, Y, 'o', ms=2.0, color='darkgray', alpha=0.7)

	ax.plot([mag_lim[0], 26.5], [0.55, 0.55], '--', color='blue', linewidth=1.5, alpha=0.75)
	ax.plot([26.5, 28.0], [0.55, -1.0], '--', color='blue', linewidth=1.5, alpha=0.75)

	cidx_cnd = ((Y < 0.) | ((Y < 0.55) & (Y < (-1.55/1.5)*(X-28.0)-1.0)))    # for cosmic ray candidates

	# Axis 3 : mag - FLUX_RADIUS
	ax = ax3
	ax.set_xticks(mag_tick)
	ax.set_xticklabels(mag_tick, fontsize=12.0)
	ax.set_yticks(flxr_tick)
	ax.set_yticklabels(flxr_tick, fontsize=12.0)
	ax.set_xlabel(band[i], fontsize=12.0)
	ax.set_ylabel('FLUX_RADIUS [pixel]', fontsize=12.0)
	ax.set_xlim(mag_lim)
	ax.set_ylim(flxr_lim)
	ax.tick_params(width=1.5, length=8.0)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
	ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
	ax.tick_params(width=1.5,length=5.0,which='minor')
	for axis in ['top','bottom','left','right']:
	    ax.spines[axis].set_linewidth(1.5)

	X = ds['mag_auto'].values
	Y = ds['flxrad'].values
	ax.plot(X, Y, 'o', ms=2.0, color='darkgray', alpha=0.7)

	ax.plot([flxr_lim[0], 26.5], [2.25, 2.25], '--', color='blue', linewidth=1.5, alpha=0.75)
	ax.plot([26.5, 28.0], [2.25, -2.0], '--', color='blue', linewidth=1.5, alpha=0.75)

	flxr_cnd = ((Y < 0.) | ((Y < 2.25) & (Y < -(4.25/1.5)*(X-28.0)-2.0)))

	# Axis 4 : mag - FWHM
	ax = ax4
	ax.set_xticks(mag_tick)
	ax.set_xticklabels(mag_tick, fontsize=12.0)
	ax.set_yticks(fwhm_tick)
	ax.set_yticklabels(fwhm_tick, fontsize=12.0)
	ax.set_xlabel(band[i], fontsize=12.0)
	ax.set_ylabel('FWHM [pixel]', fontsize=12.0)
	ax.set_xlim(mag_lim)
	ax.set_ylim(fwhm_lim)
	ax.tick_params(width=1.5, length=8.0)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
	ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
	ax.tick_params(width=1.5,length=5.0,which='minor')
	for axis in ['top','bottom','left','right']:
	    ax.spines[axis].set_linewidth(1.5)

	X = ds['mag_auto'].values
	Y = ds['fwhm'].values
	ax.plot(X, Y, 'o', ms=2.0, color='darkgray', alpha=0.7)

	ax.plot([fwhm_lim[0], 26.5], [3.0, 3.0], '--', color='blue', linewidth=1.5, alpha=0.75)
	ax.plot([26.5, 28.0], [3.0, -4.0], '--', color='blue', linewidth=1.5, alpha=0.75)

	fwhm_cnd = ((Y < 0.0) | ((Y < 3.0) & (Y < -(7.0/1.5)*(X-28.0)-4.0)))

	# Axis 5 : mag - mu0
	ax = ax5
	ax.set_xticks(mag_tick)
	ax.set_xticklabels([])
	ax.set_yticks(mu0_tick)
	ax.set_yticklabels(mu0_tick, fontsize=12.0)
	ax.set_ylabel(r'$\mu_{0}~{\rm[mag~arcsec^{-2}]}$', fontsize=12.0)
	ax.set_xlim(mag_lim)
	ax.set_ylim(mu0_lim)
	ax.tick_params(width=1.5, length=8.0)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
	ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
	ax.tick_params(width=1.5,length=5.0,which='minor')
	for axis in ['top','bottom','left','right']:
	    ax.spines[axis].set_linewidth(1.5)

	X = ds['mag_auto'].values
	Y = ds['mu0'].values
	ax.plot(X, Y, 'o', ms=2.0, color='darkgray', alpha=0.7)	

	ax.plot(mag_lim, np.array(mag_lim)-4.0, '--', color='blue', linewidth=1.5, alpha=0.75)

	mu0_cnd = ((Y < 28.5) & (Y < X-4.0))

	# Axis 6 : color-magnitude diagram
	ax = ax6
	ax.set_xticks(mag_tick)
	ax.set_xticklabels(mag_tick, fontsize=12.0)
	ax.set_yticks(col_tick)
	ax.set_yticklabels(col_tick, fontsize=12.0)
	ax.set_xlabel(band[i], fontsize=12.0)
	ax.set_xlim(mag_lim)
	ax.set_ylim(col_lim)
	ax.tick_params(width=1.5, length=8.0)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
	ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
	ax.tick_params(width=1.5,length=5.0,which='minor')
	for axis in ['top','bottom','left','right']:
	    ax.spines[axis].set_linewidth(1.5)

	X = ds['mag_auto'].values
	if (i < np.arange(len(band))[-1]):
		exec("col = ds['mag_auto']-d_sep1_"+band[i+1][1:4]+"['mag_auto']")
		ax.set_ylabel(band[i]+r'$-$'+band[i+1], fontsize=12.0)
	elif (i == np.arange(len(band))[-1]):
		exec("col = d_sep1_"+band[i-1][1:4]+"['mag_auto']-ds['mag_auto']")
		ax.set_ylabel(band[i-1]+r'$-$'+band[i], fontsize=12.0)
	Y = col.values
	ax.plot(X, Y, 'o', ms=2.0, color='darkgray', alpha=0.7)		

	col_cnd = ((Y < -0.5) | (Y > 2.5))

	plt.savefig('fig_sep1_'+band[i]+'.pdf')
	plt.savefig('fig_sep1_'+band[i]+'.png', dpi=300)
	plt.close()


# ----- Selecting cosmic ray candidates ----- #
	cos1_cnd = (mag_cnd | col_cnd | mu0_cnd | (1*cidx_cnd+1*flxr_cnd+1*fwhm_cnd >= 2))
	exec("cos1_"+band[i][1:4]+" = cos1_cnd")

	f = open("Regions/cos1_"+band[i][1:4]+".reg","w")
	for i in np.arange(np.sum(cos1_cnd)):
		f.write('{0:.2f}  {1:.2f}  {2:d}\n'.format(ds.loc[cos1_cnd]['x'].values[i],
			                                       ds.loc[cos1_cnd]['y'].values[i],
			                                       ds.loc[cos1_cnd]['num'].values[i]))
	f.close()

cos2 = (cos1_435 | cos1_606 | cos1_814)
f = open("Regions/cos2.reg","w")
for i in np.arange(np.sum(cos2)):
	f.write('{0:.2f}  {1:.2f}  {2:d}\n'.format(ds.loc[cos2]['x'].values[i],
		                                       ds.loc[cos2]['y'].values[i],
		                                       ds.loc[cos2]['num'].values[i]))
f.close()



# Printing the running time  
print("--- %s seconds ---" % (time.time() - start_time))


