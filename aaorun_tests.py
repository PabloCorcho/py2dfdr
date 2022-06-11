from py2dfdr import aaorun_command
import os
from astropy.io import fits
import numpy as np

data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', '20180310', 'ccd_1')
print(data_dir)
###################
# Reduce dark files
###################
# dark_id = ['119', '120', '121']
# dark_data = []
# for dark in dark_id:
#     path_to_dark = os.path.join(data_dir, '10mar10{}.fits'.format(dark))
#     print('{} exists?: '.format(path_to_dark), os.path.isfile(path_to_dark))
#     aaorun_command('reduce_dark', path_to_dark, idx_file='koala_dark.idx', output='dark_reduction_log.txt')
#     dark_data.append(fits.getdata(path_to_dark.replace('.fits', 'red.fits')))
# # Combine all data
# darks_to_combine = ''
# for dark in dark_id:
#     path_to_dark = os.path.join(data_dir, '10mar10{}red.fits'.format(dark))
#     darks_to_combine = ' '.join((darks_to_combine, path_to_dark))
# darks_to_combine = darks_to_combine[1:]
# # create master dark
# masterdark_name = os.path.join(os.path.dirname(path_to_dark), 'masterdark.fits')
# print('Combining darks into single master dark file: {}'.format(masterdark_name))
# aaorun_command(command='combine_image',
#                idx_file='koala_dark.idx',
#                file='\"' + darks_to_combine + '\"',
#                options=['-COMBINEDFILE %s' % masterdark_name],
#                wdir=data_dir)

############################
# Long-slit / detector  flat
############################
# lflat_id = ['031', '032', '033']
# lflat_data = []
# for lflat in lflat_id:
#     path_to_lflat = os.path.join(data_dir, '10mar10{}.fits'.format(lflat))
#     print('{} exists?: '.format(path_to_lflat), os.path.isfile(path_to_lflat))
#     aaorun_command('reduce_lflat', path_to_lflat, idx_file='koala_dark.idx', output='lflat_reduction_log.txt')
#     lflat_data.append(fits.getdata(path_to_lflat.replace('.fits', 'red.fits')))
# # Combine all data
# lflats_to_combine = ''
# for lflat in lflat_id:
#     path_to_lflat = os.path.join(data_dir, '10mar10{}red.fits'.format(lflat))
#     lflats_to_combine = ' '.join((lflats_to_combine, path_to_lflat))
# lflats_to_combine = lflats_to_combine[1:]
# # create master lflat
# masterlflat_name = os.path.join(os.path.dirname(path_to_lflat), 'masterlflat.fits')
# print('Combining lflats into single master lflat file: {}'.format(masterlflat_name))
# aaorun_command(command='combine_image',
#                idx_file='koala_dark.idx',
#                file='\"' + lflats_to_combine + '\"',
#                options=['-COMBINEDFILE %s' % masterlflat_name],
#                wdir=data_dir)

###################
# Create tram maps
###################
fibreflat = '061'
path_to_fibreflat = os.path.join(data_dir, '10mar10{}.fits'.format(fibreflat))
print('{} exists?: '.format(path_to_fibreflat), os.path.isfile(path_to_fibreflat))
aaorun_command('make_tlm', path_to_fibreflat, idx_file='koala.idx',
               output='tram_reduction_log.txt')
