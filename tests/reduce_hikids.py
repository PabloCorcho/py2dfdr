#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 07:59:41 2022

@author: pablo
"""

import yaml
import os
from py2dfdr import aaorun_command
import reduction_QC as QC
import koala_cheatsheet as kcs
import numpy as np

# =============================================================================
# PARAMS
# =============================================================================
obs_run_path = os.path.join('/media', 'pablo', 'toshiba-pab',
                            'reduce_koala_april',
                            'obs_run_0')
ccds = ['ccd_1', 'ccd_2']
dark_idx_file = 'koala_dark.idx'
fflat_idx_file = 'koala_dark.idx'
fibreflat_idx_file = 'koala_dark.idx'
arcs_idx_file = 'koala_arcs.idx'
# =============================================================================
# LOG FILE
# =============================================================================
log = open(os.path.join(obs_run_path, 'reduction_log.txt'), 'w')
log.write('Starting reduction process... \n\n')
# =============================================================================
# Load obs run info
# =============================================================================
with open(os.path.join(obs_run_path, "obs_run_info.yml"), "r") as stream:
    try:
        obs_run_info = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


nights = list(obs_run_info.keys())
[nights.remove(i) for i in ['bias', 'darks', 'fflats']]

# =============================================================================
# Reduce bias
# =============================================================================
log.write('\nreducing bias...')
for ccd in ccds:
    for name, bias_file in obs_run_info['bias'][ccd].items():
        log.write('BIAS: %s\n' % bias_file['PATH'])
        print(bias_file)
# =============================================================================
# Reduce darks
# =============================================================================
# log.write('\nStarting dark files reduction...\n')
# for ccd in ccds:
#     for exptime in obs_run_info['darks'][ccd].keys():
#         all_darks = []
#         for name, dark_file in obs_run_info['darks'][ccd][exptime].items():
#             print('·[{}] [{}] DARK: {}\n'.format(ccd, exptime,
#                                                  dark_file['PATH']))
#             log.write('·[{}] [{}] DARK: {}\n'.format(ccd, exptime,
#                                                      dark_file['PATH']))
#             path_to_dark = os.path.join(obs_run_path, 'darks', ccd, exptime,
#                                         dark_file['PATH'])
#             aaorun_command('reduce_dark', path_to_dark, idx_file=dark_idx_file,
#                            output=path_to_dark.replace('.fits', '_log.txt'))
#             all_darks.append(path_to_dark.replace('.fits', 'red.fits'))
#         print('# Computing MASTERDARK [{}] [{}]\n'.format(ccd, exptime))
#         log.write('# Computing MASTERDARK [{}] [{}]\n'.format(ccd, exptime))
#         darks_to_combine = ' '.join(all_darks)
#         masterdark_name = os.path.join(obs_run_path, 'darks', ccd, exptime,
#                                        'DARKcombined_%s.fits' % exptime)
#         logmaster = masterdark_name.replace('.fits', '_log.txt')
#         aaorun_command(command='combine_image',
#                        file='\"' + darks_to_combine + '\"',
#                        options=['-COMBINEDFILE %s' % masterdark_name],
#                        output=logmaster)
#         print('---> MASTERDARK file saved as %s \n' % masterdark_name)
#         log.write('---> MASTERDARK file saved as %s \n' % masterdark_name)
#         print('---> MASTERDARK log saved as %s \n' % logmaster)
#         log.write('---> MASTERDARK log saved as %s \n' % logmaster)
#         mdark_fig, mdark_pcnt = QC.check_image(masterdark_name)
#         mdark_fig.savefig(os.path.join(obs_run_path, 'darks', ccd, exptime,
#                                        'masterdark.png'), bbox_inches='tight')
#         log.write('---> QC plot saved as masterdark.png \n')

master_darks = {
    'ccd_1': os.path.join(obs_run_path, 'darks', 'ccd_1', '1800.0',
                          'DARKcombined_1800.0.fits'),
    'ccd_2': os.path.join(obs_run_path, 'darks', 'ccd_2', '1800.0',
                          'DARKcombined_1800.0.fits')}

# =============================================================================
# Reduce long-slit (a.k.a. detector or floppy) flats
# =============================================================================
log.write('\nStarting long-slit flat files reduction...\n')
master_lflats = {}
for ccd in ccds:
    gratings = obs_run_info['fflats'][ccd].keys()
    master_lflats[ccd] = {}
    for grating in gratings:
        master_lflats[ccd][grating] = {}
        for name, fflat_file in obs_run_info['fflats'][ccd][grating].items():
            fflat_exptime = fflat_file['EXPTIME']
            if fflat_exptime not in master_lflats[ccd][grating].keys():
                master_lflats[ccd][grating][fflat_exptime] = []
            print('·[{}] [{}] LFLAT: {}\n'.format(ccd, grating,
                                                  fflat_file['PATH']))
            log.write('·[{}] [{}] LFLAT: {}\n'.format(ccd, grating,
                                                      fflat_file['PATH']))
            path_to_fflat = os.path.join(obs_run_path, 'fflats', ccd, grating,
                                         fflat_file['PATH'])
            aaorun_command('reduce_lflat', path_to_fflat,
                           idx_file=fflat_idx_file,
                           output=path_to_fflat.replace('.fits', '_log.txt'))
            master_lflats[ccd][grating][fflat_exptime].append(
                path_to_fflat.replace('.fits', 'red.fits'))
        for exptime in master_lflats[ccd][grating].keys():
            print('# Computing MASTERLFLAT [{}] [{}] [{}]\n'.format(
                ccd, grating, exptime))
            log.write('# Computing MASTERLFLAT [{}] [{}] [{}]\n'.format(
                ccd, grating, exptime))
            fflats_to_combine = ' '.join(master_lflats[ccd][grating][exptime])
            masterfflat_name = os.path.join(
                obs_run_path, 'fflats', ccd, grating,
                'LFLATcombined_{}_{}.fits'.format(grating, exptime))
            logmaster = masterfflat_name.replace('.fits', '_log.txt')
            aaorun_command(command='combine_image',
                           file='\"' + fflats_to_combine + '\"',
                           options=['-COMBINEDFILE %s' % masterfflat_name],
                           output=logmaster)
            master_lflats[ccd][grating][exptime] = masterfflat_name
            print('---> MASTERLFLAT file saved as %s \n' % masterfflat_name)
            log.write('---> MASTERLFLAT file saved as %s \n' % masterfflat_name
                      )
            print('---> MASTERLFLAT log saved as %s \n' % logmaster)
            log.write('---> MASTERLFLAT log saved as %s \n' % logmaster)
            mdark_fig, mdark_pcnt = QC.check_image(masterfflat_name)
            mdark_fig.savefig(os.path.join(
                obs_run_path, 'fflats', ccd, grating,
                'masterlflat_{}.png'.format(exptime)), bbox_inches='tight')
            log.write('---> QC plot saved as masterlflat.png \n')
        recommended_exp_time = kcs.lflat_time[ccd][grating.split('_')[0]]
        time_keys = list(master_lflats[ccd][grating].keys())
        times = np.array(time_keys, dtype=float)
        best = np.argmin(np.abs(times - recommended_exp_time))
        master_lflats[ccd][grating] = master_lflats[ccd][grating][
            time_keys[best]]
        recommended = os.path.join(
            os.path.dirname(master_lflats[ccd][grating]), 'RECOMMENDED_FLAT')
        with open(recommended, 'w') as f:
            f.write(master_lflats[ccd][grating])


# # =============================================================================
# # Extrac tramlines
# # =============================================================================
log.write('-'*50 + '\nStarting tramline extraction from fibre flats...\n'
           + '-'*50 + '\n')
tlm_maps = {}
for night in nights:
    tlm_maps[night] = {}
    for ccd in ccds:
        tlm_maps[night][ccd] = {}
        gratings = obs_run_info[night][ccd].keys()
        for grating in gratings:
            exptimes = []
            names = []
            for name, tram_file in obs_run_info[night][ccd][grating][
                    'fibreflat'].items():
                print('· [{}] [{}] [{}]'.format(night, ccd, grating)
                      + ' TRAM: %s\n' % tram_file['PATH'])
                log.write('· [{}] [{}] [{}]'.format(night, ccd, grating)
                          + ' TRAM: %s\n' % tram_file['PATH'])
                path_to_fibreflat = os.path.join(obs_run_path, night, ccd,
                                                  grating, 'fibreflat',
                                                  tram_file['PATH'])
                aaorun_command('make_tlm', path_to_fibreflat,
                    # options=['-DARK_FILENAME %s' % master_darks[ccd]],
                    idx_file=fibreflat_idx_file,
                    output=path_to_fibreflat.replace('.fits', '_log.txt')
                    )
                exptimes.append(float(tram_file['EXPTIME']))
                names.append(name)
            recommended_exp_time = kcs.lflat_time[ccd][grating.split('_')[0]]
            best = np.argmin(np.abs(np.array(exptimes) - recommended_exp_time))
            best_fibreflat = os.path.join(
                obs_run_path, night, ccd,
                grating, 'fibreflat',
                obs_run_info[night][ccd][grating]['fibreflat'][
                    names[best]]['PATH'].replace('.fits', 'tlm.fits'))

            print('· [{}] [{}] [{}]'.format(night, ccd, grating)
                  + ' BEST RECOMMENDED TRAM: %s\n' % best_fibreflat)
            log.write('· [{}] [{}] [{}]'.format(night, ccd, grating)
                      + ' BEST RECOMMENDED TRAM: %s\n' % best_fibreflat)
            tlm_maps[night][ccd][grating] = best_fibreflat

# # # =============================================================================
# # # Reduce arcs
# # # =============================================================================
log.write('-'*50 + '\nReducing calibration arcs...\n'
          + '-'*50 + '\n')
arcs_selected = {}
for night in nights:
    for ccd in ccds:
        arcs_selected[night][ccd] = {}
        gratings = obs_run_info[night][ccd].keys()
        for grating in gratings:
            for name, arc_file in obs_run_info[night][ccd][grating][
                    'arcs'].items():
                print('· [{}] [{}] [{}] [{}]'.format(night, ccd, grating,
                                                     arc_file['ARCNAME'])
                      + ' ARC: %s\n' % arc_file['PATH'])
                log.write('· [{}] [{}] [{}] [{}]'.format(night, ccd, grating,
                                                         arc_file['ARCNAME'])
                          + ' ARC: %s\n' % arc_file['PATH'])

                path_to_arc = os.path.join(obs_run_path, night, ccd,
                                           grating, 'arcs',
                                           arc_file['PATH'])
                aaorun_command('reduce_arc', path_to_arc,
                               idx_file=arcs_idx_file,
                               options=['-DARK_FILENAME {}'.format(
                                        master_darks[ccd]),
                                        '-TLMAP_FILENAME {}'.format(
                                        tlm_maps[night][ccd][grating])],
                               output=path_to_arc.replace('.fits', '_log.txt')
                               )
                # reduced file name
                all_darks.append(os.path.join(
                    obs_run_path, 'darks',
                    ccd, exptime, dark_file['PATH'].replace('.fits', 'red.fits')))
        darks_to_combine = ' '.join(all_darks[:3])
        masterdark_name = os.path.join(obs_run_path, 'darks', ccd, exptime,
                                        'DARKcombined_%s.fits' % exptime)
        aaorun_command(command='combine_image',
                        file='\"' + darks_to_combine + '\"',
                        # options=['-COMBINEDFILE %s' % masterdark_name]
                        )
        print(darks_to_combine)


# # # Mr Krtxo
