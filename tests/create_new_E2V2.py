"""
This scripts creates a modified version of the bad pixel mask by default provided
in 2dfdr
"""
from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt

path_to_2dfdr = '/home/pablo/Software/2dfdr_install'
data = path_to_2dfdr + '/share/2dfdr/E2V2A.fits'

hdul = fits.open(data)
# Save the old version
old_data = hdul[0].data.copy()
#hdul.writeto(path_to_2dfdr + '/share/2dfdr/E2V2A_old.fits')
old_data = fits.getdata(path_to_2dfdr + '/share/2dfdr/E2V2A_old.fits')



# Create a new version
hdul[0].data = np.zeros_like(old_data)
hdul.writeto(data, overwrite=True)

plt.figure(figsize=(15, 5))
plt.subplot(121)
plt.imshow(old_data, vmax=1, vmin=0, cmap='Accent', interpolation='none')
plt.colorbar()
plt.subplot(122)
plt.imshow(hdul[0].data, vmax=1, vmin=0, cmap='Accent', interpolation='none')
plt.colorbar()
plt.show()

