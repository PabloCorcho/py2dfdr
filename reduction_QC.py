#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 08:39:39 2022

@author: pablo
"""

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt

def check_image(path, percentiles=[1, 5, 16, 50, 84, 95, 99]):
    """blah."""
    master = fits.getdata(path)
    percents = np.nanpercentile(master.flatten(), percentiles)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(121)
    if not np.isfinite(percents).all():
        print('WARNING: ALL PIXELS HAVE NAN COUNTS')
        return fig, percents
    else:
        ax.hist(master.flatten(), range=[percents[0], percents[-1]],
                bins=master.flatten().size//1000, log=True, color='k')
        ax.set_xlabel('counts')
        ax.set_ylabel('# pixels')
        for pcnt in percents:
            ax.axvline(pcnt, ls='-', color='r')
        ax = fig.add_subplot(122)
        mappable = ax.imshow(master, vmin=percents[0], vmax=percents[-1],
                             cmap='nipy_spectral', aspect='auto',
                             origin='lower')
        plt.colorbar(mappable, label='counts')
        return fig, percents


def check_tramline(path, plot=True):
    """blah."""
    with fits.open(path) as f:
        median_fwhm = f[0].header['MWIDTH']
        fibrepos = f[0].data
        mean_amplitude = f[1].data
    # Dispersion in pixels along the spectral axis
    fibrepos_disp = np.std(fibrepos, axis=1)
    intensity_disp = np.std(mean_amplitude, axis=0)
    # Quality assesment
    fibre_separation = (median_fwhm > 2) & (median_fwhm < 3)
    smoothness = ((fibrepos_disp.max() < 10) & (fibrepos_disp.min() > 0)
                  & (intensity_disp.max() < 10) & (intensity_disp.min() > 0))
    if fibre_separation & smoothness:
        bad_tramline = False
    else:
        print('[QC] · Bad tramline')
        bad_tramline = True
    # Quality control plots
    if plot:
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(6, 6),
                                gridspec_kw=dict(wspace=0.5, hspace=0.5))
        fig.suptitle(path, fontsize=8)
        ax = axs[0, 0]
        ax.set_title(r'$\mu(\lambda)$')
        mappable = ax.imshow(fibrepos, aspect='auto', origin='lower')
        plt.colorbar(mappable, ax=ax)
        ax.set_ylabel('Fibre')
        ax.set_xlabel(r'$\lambda (pix)$')
        ax = axs[0, 1]
        ax.set_title('Median FWHM\nbetween fibres: {:.2f} pix'.format(
            median_fwhm))
        ax.plot(fibrepos_disp)
        ax.set_ylabel(r'$\sigma(\mu)$')
        ax.set_xlabel('Fibre')
        ax = axs[1, 0]
        ax.set_title(r'$I(\lambda)$')
        mappable = ax.imshow(mean_amplitude, aspect='auto', origin='lower')
        plt.colorbar(mappable, ax=ax)
        ax.set_ylabel('Fibre')
        ax.set_xlabel(r'$\lambda (pix)$')
        ax = axs[1, 1]
        ax.plot(intensity_disp)
        ax.set_ylabel(r'$\sigma(I)$')
        ax.set_xlabel(r'$\lambda (pix)$')
        fig.savefig(path.replace('.fits', '.png'), bbox_inches='tight')
        plt.clf()
        plt.close(fig)
    return bad_tramline


def clean_nan(path):
    """blah..."""
    with fits.open(path, mode='update') as hdul:
        nans_mask = np.isfinite(hdul[0].data)
        fits_shape = nans_mask.shape
        nans_pos = np.where(nans_mask == False)
        if nans_pos[0].size > 0:
            for entry in zip(nans_pos[0], nans_pos[1]):
                hdul[0].data[entry[0], entry[1]] = np.nanmean([
                    hdul[0].data[np.clip(entry[0]-1, a_min=0, a_max=None),
                                 entry[1]],
                    hdul[0].data[entry[0],
                                 np.clip(entry[1]-1, a_min=0, a_max=None)],
                    hdul[0].data[np.clip(entry[0]+1, a_min=None,
                                         a_max=fits_shape[0]-1),
                                 entry[1]],
                    hdul[0].data[entry[0],
                                 np.clip(entry[1]+1, a_min=None,
                                         a_max=fits_shape[1]-1)]])
            hdul.flush()
        hdul.close()

# Mr. Krtxo
