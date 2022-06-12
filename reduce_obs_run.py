#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 15:50:34 2022

@author: pablo
"""

import yaml
import os
from py2dfdr import aaorun_command
import reduction_QC as QC
import koala_cheatsheet as kcs
import numpy as np
from matplotlib import pyplot as plt
import logging
import datetime


class ReduceObsRun(object):
    """Reduce KOALA observing runs."""

    def __init__(self, obs_run_path, **kwargs):
        self.obs_run_path = obs_run_path
        # LOG FILE ------------------------------------------------------------
        logging.basicConfig(
            filename=os.path.join(self.obs_run_path, 'OR_reduction.log'),
            level=logging.INFO)
        logging.info(
            '-'*50 + '\n'
            'Initialising OR reduction process at:\n   {} \n'.format(
                self.obs_run_path) + '-'*50 + '\n')
        logging.info(datetime.datetime.now().strftime("%c"))
        # CCD detector to reduce ----------------------------------------------
        self.ccds = kwargs.get('ccds', None)
        logging.info('CCDs data to reduce: {} \n'.format(', '.join(self.ccds)))
        if self.ccds is None:
            logging.error('CCDs not provided [ccd_1, ccd_2] \n')
            raise NameError('CCDs not provided [ccd_1, ccd_2]')

        # Parameter files for 2DFDR -------------------------------------------
        logging.info('-> Setting configuration files for 2dfdr \n')
        self.dark_idx_file = kwargs.get('dark_idx', None)
        logging.info('--> DARK configuration file: %s\n' % self.dark_idx_file)
        self.lflat_idx_file = kwargs.get('lflat_idx', None)
        logging.info('--> LFLAT configuration file: %s\n' % self.lflat_idx_file)
        self.fibreflat_idx_file = kwargs.get('fibreflat_idx', None)
        logging.info('--> FFLAT configuration file: %s\n' % self.fibreflat_idx_file)
        self.arc_idx_file = kwargs.get('arcs_idx', None)
        logging.info('--> ARCS configuration file: %s\n' % self.arc_idx_file)
        self.object_idx_file = kwargs.get('object_idx', None)
        logging.info('--> OBJECT configuration file: %s\n' % self.object_idx_file)

        # Load yml file containing data description ---------------------------
        self.load_obs_run_yml()

    def load_obs_run_yml(self):
        """Load yaml file of the Observing Run."""
        logging.info('Loading yml file containing the OR data\n')
        with open(os.path.join(self.obs_run_path, "obs_run_info.yml"),
                  "r") as stream:
            try:
                self.obs_run_info = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                logging.error('· [ERROR] Unable to load yml file\n')
                logging.error(exc)
        # Get all the nights contained within the observing run
        self.nights = list(self.obs_run_info.keys())
        [self.nights.remove(i) for i in ['bias', 'darks', 'fflats']]

    def reduce_bias(self):
        """blah."""
        logging.info('\nreducing bias...')
        for ccd in self.ccds:
            for name, bias_file in self.obs_run_info['bias'][ccd].items():
                logging.info('BIAS: %s\n' % bias_file['PATH'])

    def reduce_darks(self):
        """blah."""
        logging.info('-'*50 + '\n    Reducing DARKS    \n' + '-'*50 + '\n')
        for ccd in self.ccds:
            for exptime in self.obs_run_info['darks'][ccd].keys():
                all_darks = []
                for name, dark_file in self.obs_run_info[
                        'darks'][ccd][exptime].items():
                    logging.info('·[{}] [{}] DARK: {}\n'.format(
                        ccd, exptime,
                        dark_file['PATH']))

                    path_to_dark = os.path.join(
                        self.obs_run_path, 'darks', ccd, exptime,
                        dark_file['PATH'])
                    aaorun_command('reduce_dark', path_to_dark,
                                   idx_file=self.dark_idx_file,
                                   output=path_to_dark.replace(
                                       '.fits', '_log.txt'))
                    all_darks.append(path_to_dark.replace('.fits', 'red.fits'))
                logging.info(
                    '# Computing MASTERDARK [{}] [{}]\n'.format(ccd, exptime))
                darks_to_combine = ' '.join(all_darks)
                masterdark_name = os.path.join(
                    self.obs_run_path, 'darks',
                    ccd, exptime, 'DARKcombined_%s.fits' % exptime)
                logmaster = masterdark_name.replace('.fits', '_log.txt')
                # CALL TO AAORUN
                aaorun_command(command='combine_image',
                               file='\"' + darks_to_combine + '\"',
                               options=['-COMBINEDFILE %s' % masterdark_name],
                               output=logmaster)
                logging.info(
                    '---> MASTERDARK file saved as %s \n' % masterdark_name)
                logging.info('---> MASTERDARK log saved as %s \n' % logmaster)
                # Sanity check plots
                mdark_fig, mdark_pcnt = QC.check_image(masterdark_name)
                mdark_fig.savefig(os.path.join(self.obs_run_path, 'darks',
                                               ccd, exptime, 'masterdark.png'),
                                  bbox_inches='tight')
                plt.clf()
                plt.close()
                logging.info('---> QC plot saved as masterdark.png \n')

        self.master_darks = {
            'ccd_1': os.path.join(self.obs_run_path, 'darks', 'ccd_1',
                                  '1800.0', 'DARKcombined_1800.0.fits'),
            'ccd_2': os.path.join(self.obs_run_path, 'darks', 'ccd_2',
                                  '1800.0', 'DARKcombined_1800.0.fits')}

    def get_master_darks(self, exptime='1800.0'):
        """blah."""
        self.master_darks = {}
        print(' --> Searching master dark files\n')
        logging.info(' --> Searching master dark files\n')
        for ccd in self.ccds:
            path_to_master = os.path.join(
                self.obs_run_path, 'darks', ccd, exptime,
                'DARKcombined_{}.fits'.format(exptime))
            if os.path.isfile(path_to_master):
                print('--> [{}] MASTERDARK found at {}\n'.format(
                        ccd, path_to_master))
                logging.info('--> [{}] MASTERDARK found at {}\n'.format(
                    ccd, path_to_master))
                self.master_darks[ccd] = path_to_master
            else:
                logging.error(
                    '--> [{}] [ERROR] MASTERDARK *NOT* found at {}\n'.format(
                        ccd, path_to_master))
                raise NameError('[ERROR] Unable to find a master dark (see log)')

    def reduce_lflats(self):
        """blah."""
        logging.info('\nStarting long-slit flat files reduction...\n')
        self.master_lflats = {}
        for ccd in self.ccds:
            gratings = self.obs_run_info['fflats'][ccd].keys()
            self.master_lflats[ccd] = {}
            for grating in gratings:
                self.master_lflats[ccd][grating] = {}
                for name, fflat_file in self.obs_run_info['fflats'][ccd][
                        grating].items():
                    fflat_exptime = fflat_file['EXPTIME']
                    if fflat_exptime not in self.master_lflats[ccd][
                            grating].keys():
                        self.master_lflats[ccd][grating][fflat_exptime] = []
                    logging.info('·[{}] [{}] LFLAT: {}\n'.format(
                        ccd, grating, fflat_file['PATH']))
                    path_to_fflat = os.path.join(
                        self.obs_run_path, 'fflats', ccd, grating,
                        fflat_file['PATH'])
                    # CALL TO AAORUN
                    aaorun_command('reduce_lflat', path_to_fflat,
                                   idx_file=self.lflat_idx_file,
                                   output=path_to_fflat.replace(
                                       '.fits', '_log.txt'))
                    self.master_lflats[ccd][grating][fflat_exptime].append(
                        path_to_fflat.replace('.fits', 'red.fits'))
                for exptime in self.master_lflats[ccd][grating].keys():
                    logging.info(
                        '# Computing MASTERLFLAT [{}] [{}] [{}]\n'.format(
                            ccd, grating, exptime))
                    fflats_to_combine = ' '.join(
                        self.master_lflats[ccd][grating][exptime])
                    masterfflat_name = os.path.join(
                        self.obs_run_path, 'fflats', ccd, grating,
                        'LFLATcombined_{}_{}.fits'.format(grating, exptime))
                    logmaster = masterfflat_name.replace('.fits', '_log.txt')
                    # CALL TO AAORUN
                    aaorun_command(command='combine_image',
                                   file='\"' + fflats_to_combine + '\"',
                                   options=['-COMBINEDFILE %s' % masterfflat_name],
                                   output=logmaster)
                    self.master_lflats[ccd][grating][exptime] = masterfflat_name

                    logging.info(
                        '---> MASTERLFLAT file saved as %s \n' % masterfflat_name
                              )
                    logging.info('---> MASTERLFLAT log saved as %s \n' % logmaster)
                    mfflat_fig, mfflat_pcnt = QC.check_image(masterfflat_name)
                    mfflat_fig.savefig(os.path.join(
                        self.obs_run_path, 'fflats', ccd, grating,
                        'masterlflat_{}.png'.format(exptime)), bbox_inches='tight')
                    plt.clf()
                    plt.close()
                    logging.info('---> QC plot saved as masterlflat.png \n')
                recommended_exp_time = kcs.lflat_time[ccd][grating.split('_')[0]]
                time_keys = list(self.master_lflats[ccd][grating].keys())
                times = np.array(time_keys, dtype=float)
                best = np.argmin(np.abs(times - recommended_exp_time))
                self.master_lflats[ccd][grating] = self.master_lflats[ccd][grating][
                    time_keys[best]]
                recommended = os.path.join(
                    os.path.dirname(self.master_lflats[ccd][grating]),
                    'RECOMMENDED_FLAT')
                with open(recommended, 'w') as f:
                    f.write(self.master_lflats[ccd][grating])

    def get_master_lflats(self):
        """blah."""
        self.master_lflats = {}
        logging.info(' --> Searching master dark files\n')
        for ccd in self.ccds:
            gratings = self.obs_run_info['fflats'][ccd].keys()
            self.master_lflats[ccd] = {}
            for grating in gratings:
                recommended_file = os.path.join(
                    self.obs_run_path, 'fflats', ccd, grating,
                    'RECOMMENDED_FLAT')
                with open(recommended_file, 'r') as f:
                    path_to_master = f.readline()
                if os.path.isfile(path_to_master):
                    logging.info(
                        '--> [{}] [{}] MASTERLFLAT found at {}\n'.format(
                            ccd, grating, path_to_master))
                    self.master_lflats[ccd][grating] = path_to_master
                else:
                    logging.error(
                        '--> [{}] [{}] [ERROR] MASTERFLAT *NOT* found at {}\n'.format(
                            ccd, grating, path_to_master))
                    raise NameError('[ERROR] Unable to find a master fflat (see log)')

    def extract_tramlines(self):
        """blah."""
        logging.info('-'*50
                     + '\nStarting tramline extraction from fibre flats...\n'
                     + '-'*50 + '\n')
        self.tlm_maps = {}
        if self.fibreflat_idx_file is None:
            logging.error('--> [ERROR] No FIBREFLAT IDX file provided')
            raise NameError('--> [ERROR] No FIBREFLAT IDX file provided')
        for night in self.nights:
            self.tlm_maps[night] = {}
            for ccd in self.ccds:
                self.tlm_maps[night][ccd] = {}
                gratings = self.obs_run_info[night][ccd].keys()
                for grating in gratings:
                    exptimes = []
                    names = []
                    for name, tram_file in self.obs_run_info[night][ccd][
                            grating]['fibreflat'].items():
                        logging.info('· [{}] [{}] [{}]'.format(night, ccd,
                                                               grating)
                                     + ' TRAM: %s\n' % tram_file['PATH'])
                        path_to_fibreflat = os.path.join(
                            self.obs_run_path, night, ccd, grating,
                            'fibreflat', tram_file['PATH'])
                        # CALL TO AAORUN
                        aaorun_command('make_tlm', path_to_fibreflat,
                                       idx_file=self.fibreflat_idx_file,
                                       output=path_to_fibreflat.replace(
                                           '.fits', '_log.txt'))
                        exptimes.append(float(tram_file['EXPTIME']))
                        names.append(name)
                    recommended_exp_time = kcs.lflat_time[ccd][
                        grating.split('_')[0]]
                    best = np.argmin(
                        np.abs(np.array(exptimes) - recommended_exp_time))
                    best_fibreflat = os.path.join(
                        self.obs_run_path, night, ccd,
                        grating, 'fibreflat',
                        self.obs_run_info[night][ccd][grating]['fibreflat'][
                            names[best]]['PATH'].replace('.fits', 'tlm.fits'))

                    logging.info('· [{}] [{}] [{}]'.format(night, ccd, grating)
                                 + ' BEST RECOMMENDED TRAM: %s\n'
                                 % best_fibreflat)
                    self.tlm_maps[night][ccd][grating] = best_fibreflat
                    recommended = os.path.join(
                        os.path.dirname(best_fibreflat), 'RECOMMENDED_TRAM')
                    with open(recommended, 'w') as f:
                        f.write(best_fibreflat)

    def get_tlm_maps(self):
        """blah."""
        self.tlm_maps = {}
        logging.info(' --> Searching TRAMLINE MAPS\n')
        for night in self.nights:
            self.tlm_maps[night] = {}
            for ccd in self.ccds:
                self.tlm_maps[night][ccd] = {}
                gratings = self.obs_run_info[night][ccd].keys()
                for grating in gratings:
                    recommended = os.path.join(
                        self.obs_run_path, night, ccd, grating, 'fibreflat',
                        'RECOMMENDED_TRAM')
                    with open(recommended, 'r') as f:
                        path_to_tram = f.readline()
                    if os.path.isfile(path_to_tram):
                        logging.info(
                            '--> [{}] [{}] [{}] TRAM MAP found at {}\n'.format(
                                night, ccd, grating, path_to_tram))
                        self.tlm_maps[night][ccd][grating] = path_to_tram
                    else:
                        logging.error(
                            '--> [{}] [{}] [{}] [ERROR] TRAM MAP *NOT* found at {}\n'.format(
                                night, ccd, grating, path_to_tram))
                        raise NameError('[ERROR] Unable to find a TRAM MAP (see log)')

    def reduce_arcs(self):
        """blah."""
        logging.info('-' * 50 + '\nReducing calibration arcs...\n' + '-' * 50
                     + '\n')
        self.arcs_selected = {}
        try:
            self.master_darks
        except Exception:
            self.get_master_darks()
        try:
            self.master_lflats
        except Exception:
            self.get_master_lflats()
        try:
            self.tlm_maps
        except Exception:
            self.get_tlm_maps()

        for night in self.nights:
            self.arcs_selected[night] = {}
            for ccd in self.ccds:
                self.arcs_selected[night][ccd] = {}
                gratings = self.obs_run_info[night][ccd].keys()
                for grating in gratings:
                    exptimes = []
                    names = []
                    arcnames = []
                    for name, arc_file in self.obs_run_info[night][ccd][
                            grating]['arcs'].items():
                        arc_name = arc_file['ARCNAME']
                        arc_exptime = arc_file['EXPTIME']
                        logging.info('· [{}] [{}] [{}] [{}] [{}]'.format(
                            night, ccd, grating, arc_name, arc_exptime)
                              + ' ARC: %s\n' % arc_file['PATH'])
                        path_to_arc = os.path.join(
                            self.obs_run_path, night, ccd, grating, 'arcs',
                            arc_file['PATH'])
                        # CALL TO AAORUN
                        aaorun_command(
                            'reduce_arc', path_to_arc,
                            idx_file=self.arc_idx_file,
                            options=['-DARK_FILENAME {}'.format(
                                self.master_darks[ccd]),
                                '-TLMAP_FILENAME {}'.format(
                                    self.tlm_maps[night][ccd][grating])],
                            output=path_to_arc.replace('.fits', '_log.txt')
                            )
                        arc_fig, arc_pcnt = QC.check_image(
                            path_to_arc.replace('.fits', 'red.fits'))
                        arc_fig.savefig(
                            os.path.join(os.path.dirname(path_to_arc),
                                         'arc_{}_{}.png'.format(arc_name,
                                                                arc_exptime)),
                            bbox_inches='tight')
                        plt.clf()
                        plt.close()

                        names.append(name)
                        exptimes.append(arc_exptime)
                        arcnames.append(arc_name)

                    rec_arc_name, rec_exptime = kcs.arc_time[ccd][
                        grating.split('_')[0]]
                    selected_arcs = np.array(
                        [rec_arc_name in name for name in arcnames], dtype=bool
                        )
                    exptimes = np.array(exptimes, dtype=float)
                    exptimes[~selected_arcs] = np.nan
                    best = np.nanargmin(
                        np.abs(exptimes[selected_arcs] - rec_exptime))
                    best_arc = os.path.join(
                        self.obs_run_path, night, ccd,
                        grating, 'arcs',
                        self.obs_run_info[night][ccd][grating]['arcs'][
                            names[best]]['PATH'].replace('.fits', 'red.fits'))
                    logging.info(
                        '· [{}] [{}] [{}]'.format(night, ccd, grating)
                        + ' BEST RECOMMENDED TRAM: %s\n' % best_arc)
                    self.arcs_selected[night][ccd][grating] = best_arc
                    recommended = os.path.join(
                        os.path.dirname(best_arc), 'RECOMMENDED_ARC')
                    with open(recommended, 'w') as f:
                        f.write(best_arc)

    def get_arcs(self):
        """blah."""
        self.arcs_selected = {}
        logging.info(' --> Searching ARCS\n')
        for night in self.nights:
            self.arcs_selected[night] = {}
            for ccd in self.ccds:
                self.arcs_selected[night][ccd] = {}
                gratings = self.obs_run_info[night][ccd].keys()
                for grating in gratings:
                    recommended = os.path.join(
                        self.obs_run_path, night, ccd, grating, 'arcs',
                        'RECOMMENDED_ARC')
                    with open(recommended, 'r') as f:
                        path_to_arc = f.read().split('\n')[0]
                    if os.path.isfile(path_to_arc):
                        logging.info(
                            '--> [{}] [{}] [{}] ARC found at {}\n'.format(
                                night, ccd, grating, path_to_arc))
                        self.arcs_selected[night][ccd][grating] = path_to_arc
                    else:
                        logging.error(
                            '--> [{}] [{}] [{}] [ERROR] ARC *NOT* found at {}\n'.format(
                                night, ccd, grating, path_to_arc))
                        raise NameError('Unable to find ARC \n{}\n (see log)'.
                                        format(path_to_arc))

    def reduce_fflats(self):
        """blah."""
        logging.info('-' * 50 + '\nReducing FIBRE FLATS...\n' + '-' * 50
                     + '\n')
        try:
            self.master_darks
        except Exception:
            self.get_master_darks()
        try:
            self.master_lflats
        except Exception:
            self.get_master_lflats()
        try:
            self.tlm_maps
        except Exception:
            self.get_tlm_maps()
        try:
            self.arcs_selected
        except Exception:
            self.get_arcs()

        self.fibreflat_selected = self.tlm_maps.copy()

        for night in self.nights:
            for ccd in self.ccds:
                gratings = self.obs_run_info[night][ccd].keys()
                for grating in gratings:
                    for name, object_file in self.obs_run_info[night][ccd][
                            grating]['fibreflat'].items():
                        exptime = object_file['EXPTIME']
                        logging.info('· [{}] [{}] [{}] [{}] '.format(
                            night, ccd, grating, exptime)
                              + ' FFLAT: %s\n' % object_file['PATH'])
                        path_to_fflat = os.path.join(self.obs_run_path, night,
                                                     ccd, grating, 'fibreflat',
                                                     object_file['PATH'])
                        # CALL TO AAORUN
                        aaorun_command(
                            'reduce_fflat', path_to_fflat,
                            idx_file=self.fibreflat_idx_file,
                            options=[
                                '-DARK_FILENAME {}'.format(
                                    self.master_darks[ccd]),
                                '-TLMAP_FILENAME {}'.format(
                                    self.tlm_maps[night][ccd][grating]),
                                '-WAVEL_FILENAME {}'.format(
                                    self.arcs_selected[night][ccd][grating])],
                            output=path_to_fflat.replace('.fits', '_log.txt')
                            )
                        fflat_fig, fflat_pcnt = QC.check_image(
                            path_to_fflat.replace('.fits', 'red.fits'))
                        fflat_fig.savefig(
                            os.path.join(os.path.dirname(path_to_fflat),
                                         '{}.png'.format(name)),
                            bbox_inches='tight')
                        plt.clf()
                        plt.close()
                    self.fibreflat_selected[night][ccd][grating] = (
                        self.fibreflat_selected[night][ccd][grating].replace(
                            'tlm.fits', 'red.fits'))

    def get_fibreflats(self):
        """blah."""
        self.fibreflat_selected = {}
        logging.info(' --> Searching FIBRE FLATS\n')
        for night in self.nights:
            self.fibreflat_selected[night] = {}
            for ccd in self.ccds:
                self.fibreflat_selected[night][ccd] = {}
                gratings = self.obs_run_info[night][ccd].keys()
                for grating in gratings:
                    recommended = os.path.join(
                        self.obs_run_path, night, ccd, grating, 'fibreflat',
                        'RECOMMENDED_TRAM')
                    with open(recommended, 'r') as f:
                        path_to_tram = f.readline()
                    path_to_fflat = path_to_tram.replace(
                        'tlm.fits', 'red.fits')
                    if os.path.isfile(path_to_fflat):
                        logging.info(
                            '--> [{}] [{}] [{}] FFLAT found at {}\n'.format(
                                night, ccd, grating, path_to_fflat))
                        self.fibreflat_selected[night][ccd][grating] = path_to_fflat
                    else:
                        logging.error(
                            '--> [{}] [{}] [{}] [ERROR] FFLAT *NOT* found at {}\n'.format(
                                night, ccd, grating, path_to_fflat))
                        raise NameError(
                            '[ERROR] Unable to find a FFLAT (see log)')

    def reduce_object(self):
        """blah."""
        logging.info('-' * 50 + '\nReducing SCIENCE OBJECTS...\n'
                       + '-' * 50 + '\n')
        try:
            self.master_darks
        except Exception:
            self.get_master_darks()
        try:
            self.master_lflats
        except Exception:
            self.get_master_lflats()
        try:
            self.tlm_maps
        except Exception:
            self.get_tlm_maps()
        try:
            self.arcs_selected
        except Exception:
            self.get_arcs()
        try:
            self.fibreflat_selected
        except Exception:
            self.get_fibreflats()

        for night in self.nights:
            for ccd in self.ccds:
                gratings = self.obs_run_info[night][ccd].keys()
                for grating in gratings:
                    for name, object_file in self.obs_run_info[night][ccd][
                            grating]['sci'].items():
                        obj_name = object_file['NAME']
                        exptime = object_file['EXPTIME']
                        if obj_name.find('FOCUS') >= 0:
                            continue
                        logging.info('· [{}] [{}] [{}] [{}] [{}]'.format(
                            night, ccd, grating, obj_name, exptime)
                              + ' OBJECT: %s\n' % object_file['PATH'])
                        path_to_obj = os.path.join(
                            self.obs_run_path, night, ccd, grating, 'sci',
                            object_file['PATH'])
                        # CALL TO AAORUN
                        aaorun_command(
                            'reduce_object', path_to_obj,
                            idx_file=self.object_idx_file,
                            options=[
                                '-DARK_FILENAME {}'.format(
                                    self.master_darks[ccd]),
                                '-FFLAT_FILENAME {}'.format(
                                    self.fibreflat_selected[night][ccd][
                                        grating]),
                                '-TLMAP_FILENAME {}'.format(
                                    self.tlm_maps[night][ccd][grating]),
                                '-WAVEL_FILENAME {}'.format(
                                    self.arcs_selected[night][ccd][grating])],
                            output=path_to_obj.replace('.fits', '_log.txt')
                            )
                        obj_fig, obj_pcnt = QC.check_image(
                            path_to_obj.replace('.fits', 'red.fits'))
                        obj_fig.suptitle(obj_name, fontsize=12)
                        obj_fig.savefig(
                            os.path.join(os.path.dirname(path_to_obj),
                                         '{}.png'.format(name)),
                            bbox_inches='tight')
                        plt.clf()
                        plt.close()

    def combine_science_flats(self):
        """Combine all dome/sky flats into a master file for each night."""
        pass


if __name__ == '__main__':

    import time
    tstart = time.time()
    obsrunpath = '/media/pablo/toshiba-pab/reduce_koala_april/obs_run_0'
    redOR = ReduceObsRun(
        obs_run_path=obsrunpath,
        ccds=['ccd_1', 'ccd_2'],
        dark_idx='koala_dark.idx',
        fibreflat_idx='koala_dark.idx',
        lflat_idx='koala_dark.idx',
        arcs_idx='koala_arcs.idx',
        object_idx='koala_reduce.idx')

    # redOR.reduce_darks()
    # redOR.reduce_lflats()
    # redOR.extract_tramlines()
    # redOR.reduce_arcs()


    redOR.get_master_darks()
    redOR.get_master_lflats()
    redOR.get_tlm_maps()
    redOR.get_arcs()
    redOR.get_fibreflats()
    # redOR.reduce_fflats()

    redOR.reduce_object()

    redOR.log.close()
    tend = time.time()
    print('\n\n ### Elapsed time (hrs): ', (tend - tstart)/3600)
