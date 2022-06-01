#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 08:39:39 2022

@author: pablo
"""

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt


def check_image(path, percentiles=[5, 16, 50, 84, 95]):
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
                bins=master.flatten().size//1000, log=True)
        ax.set_xlabel('counts')
        ax.set_ylabel('# pixels')
        for pcnt in percents:
            ax.axvline(pcnt, ls='-', color='k')
        ax = fig.add_subplot(122)
        mappable = ax.imshow(master, vmin=percents[0], vmax=percents[-1],
                             cmap='nipy_spectral', aspect='auto', origin='lower')
        plt.colorbar(mappable, label='counts')
        return fig, percents

# def check_arc(path, percentiles=[5, 16, 50, 84, 95]):
#     """blah."""
#     master = fits.getdata(path)
#     percents = np.nanpercentile(master.flatten(), percentiles)
#     fig = plt.figure(figsize=(10, 5))
#     ax = fig.add_subplot(121)
#     ax.hist(master.flatten(), range=[percents[0], percents[-1]],
#             bins=master.flatten().size//1000, log=True)
#     ax.set_xlabel('counts')
#     ax.set_ylabel('# pixels')
#     for pcnt in percents:
#         ax.axvline(pcnt, ls='-', color='k')
#     ax = fig.add_subplot(122)
#     mappable = ax.imshow(master, vmin=percents[0], vmax=percents[-1],
#                          cmap='nipy_spectral', aspect='auto', origin='lower')
#     plt.colorbar(mappable, label='counts')
#     return fig, percents

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
