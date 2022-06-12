#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 08:57:30 2022

@author: pablo
"""

from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits

hdul = fits.open(
    '/media/pablo/toshiba-pab/reduce_koala_april/obs_run_0/night_20220405/ccd_2/385R_7700/fibreflat/05apr20032tlm.fits')

hdul2 = fits.open(
    '/media/pablo/toshiba-pab/reduce_koala_april/obs_run_0/night_20220404/ccd_2/385R_7700/fibreflat/04apr20031tlm.fits')


plt.figure(figsize=(10, 5))
plt.subplot(121)
plt.imshow(hdul[0].data, aspect='auto')
plt.colorbar()
plt.subplot(122)
plt.imshow(hdul2[0].data, aspect='auto')
plt.colorbar()

plt.figure(figsize=(10, 5))
plt.subplot(121)
plt.hist(hdul[0].data.flatten(), bins='auto')
plt.subplot(122)
plt.hist(hdul2[0].data.flatten(), bins='auto')


plt.figure(figsize=(10, 5))
plt.subplot(121)
plt.imshow(hdul[1].data, aspect='auto')
plt.colorbar()
plt.subplot(122)
plt.imshow(hdul2[1].data, aspect='auto')
plt.colorbar()

plt.figure(figsize=(10, 5))
plt.subplot(121)
plt.hist(hdul[1].data.flatten(), bins='auto')
plt.subplot(122)
plt.hist(hdul2[1].data.flatten(), bins='auto')