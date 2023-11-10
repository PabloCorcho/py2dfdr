"""
This scripts creates a modified version of the bad pixel mask by default provided
in 2dfdr
"""
from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
import os

#path_to_2dfdr = os.path.dirname(os.path.dirname(os.system("which aaorun")))
#path_to_2dfdr = '/home/pablo/Software/2dfdr_install'
path_to_2dfdr = "/Users/users/caballero/DATASERVER3/KOALA/2dfdr/2dfdr_install"
print(f"2dfdr directory: {path_to_2dfdr}")

data = os.path.join(path_to_2dfdr, 'share', '2dfdr', 'E2V2A.fits')
print(f"Path to detector: {data}")

hdul = fits.open(data)
# Save the old version
old_data = hdul[0].data.copy()
#hdul.writeto(path_to_2dfdr + '/share/2dfdr/E2V2A_old.fits')
original_file = data.replace(".fits", "_old.fits")
if not os.path.isfile(original_file):
    print(f"Create a backup of the old file: {original_file}")
    hdul.writeto(original_file, overwrite=True)
old_data = fits.getdata(original_file)


# Create a new version
hdul[0].data = np.zeros_like(old_data)
hdul.writeto(data, overwrite=True)
print(f"New detector data saved as {data}")

plt.figure(figsize=(15, 5))
plt.subplot(121)
plt.imshow(old_data, vmax=1, vmin=0, cmap='Accent', interpolation='none')
plt.colorbar()
plt.subplot(122)
plt.imshow(hdul[0].data, vmax=1, vmin=0, cmap='Accent', interpolation='none')
plt.colorbar()
plt.show()

