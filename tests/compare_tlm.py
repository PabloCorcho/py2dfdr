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
    '/home/pablo/Research/obs_data/HI-KIDS/raw/mar2022_obsruns/obs_run_0/night_20220308/ccd_2/385R_7700/fibreflat/08mar20005tlm.fits')

hdul2 = fits.open(
    '/home/pablo/Research/obs_data/HI-KIDS/raw/mar2022_obsruns/obs_run_0/night_20220308/ccd_1/580V_4700/fibreflat/08mar10006tlm.fits')


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