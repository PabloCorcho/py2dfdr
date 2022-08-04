#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 15:50:34 2022

@author: pablo
"""

import yaml
import os
from .py2dfdr import aaorun_command
from . import reduction_QC as QC
from . import koala_cheatsheet as kcs
import numpy as np
import logging
import datetime
import verbose


class ReduceObsRun(object):
    """Reduce KOALA observing runs."""

    def __init__(self, obs_run_path, verb=True, **kwargs):
        """..."""
        self.master_bias = None
        self.obs_run_info = None
        self.nights = None
        self.master_darks = None
        self.master_lflats = None
        self.master_tlm = None
        self.master_arcs = None
        self.master_fibreflats = None
        # ---------------------------------------------------------------------
        self.obs_run_path = obs_run_path
        # Initialise logging file
        logging.basicConfig(
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%m/%d/%Y %I:%M:%S %p',
            filename=os.path.join(self.obs_run_path, 'OR_reduction.log'),
            # handlers=[
            #    logging.FileHandler("OR_reduction.log"),
            #    logging.StreamHandler()],
            level=logging.DEBUG)
        if verb:
            logging.getLogger().addHandler(logging.StreamHandler())
        verbose.log_header('Initialising OR reduction process at:\n   {} \n'.format(
            self.obs_run_path))
        logging.info(datetime.datetime.now().strftime("%c"))
        # CCD detector to reduce
        self.ccds = kwargs.get('ccds', None)
        if self.ccds is None:
            logging.error('CCDs not provided [ccd_1, ccd_2] \n')
            raise NameError('CCDs not provided [ccd_1, ccd_2]')
        else:
            logging.info('CCDs data to reduce: {} \n'.format(', '.join(self.ccds)))
        # Parameter files for 2dfdr
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

        # Load yml file containing data description
        self.load_obs_run_yml(**kwargs)

    def load_obs_run_yml(self, **kwargs):
        """Load the the Observing Run yaml file."""
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
        # Extra
        night_to_remove = kwargs.get('night_to_remove', [])
        [self.nights.remove(i) for i in ['bias', 'darks', 'fflats'] + night_to_remove]

    def check_masters(self, names=None):
        """Check if a dictionary containing master files is loaded in memory."""
        masters = {'bias': self.master_bias, 'darks': self.master_darks,
                   'lflats': self.master_lflats, 'tlm': self.master_tlm,
                   'arcs': self.master_arcs, 'fflats': self.master_fibreflats}
        get_method = {'bias': self.get_master_bias, 'darks': self.get_master_darks,
                      'lflats': self.get_master_lflats, 'tlm': self.get_master_tlm,
                      'arcs': self.get_master_arcs, 'fflats': self.get_master_fibreflats}
        if not names:
            pass
        else:
            for name in names:
                assert masters[name] is not None, get_method[name]()

    def reduce_bias(self):
        """blah."""
        verbose.log_header('bias')
        for ccd in self.ccds:
            for name, bias_file in self.obs_run_info['bias'][ccd].items():
                logging.info('BIAS: %s\n' % bias_file['PATH'])

    def get_master_bias(self):
        pass

    def reduce_darks(self):
        """blah."""
        verbose.log_header('Reducing darks')
        assert self.dark_idx_file is not None, verbose.missing_idx('dark_idx')
        for ccd in self.ccds:
            for exptime in self.obs_run_info['darks'][ccd].keys():
                all_darks = []
                for name, dark_file in self.obs_run_info['darks'][ccd][exptime].items():
                    logging.info('·[{}] [{}] DARK: {}\n'.format(
                        ccd, exptime,
                        dark_file['PATH']))

                    path_to_dark = os.path.join(
                        self.obs_run_path, 'darks', ccd, exptime,
                        dark_file['PATH'])
                    # CALL AAORUN
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
                # CALL AAORUN COMBINE
                aaorun_command(command='combine_image',
                               file='\"' + darks_to_combine + '\"',
                               idx_file=self.dark_idx_file,
                               options=['-COMBINEDFILE %s' % masterdark_name],
                               output=logmaster,
                               wdir=os.path.dirname(masterdark_name),
                               log=True)
                logging.info(
                    '---> MASTERDARK file saved as %s \n' % masterdark_name)
                logging.info('---> MASTERDARK log saved as %s \n' % logmaster)
                # Sanity check plots
                QC.check_image(masterdark_name,
                               save_dir=os.path.join(self.obs_run_path, 'darks',
                                                     ccd, exptime, 'masterdark.png'))
        self.master_darks = {}
        for ccd in self.ccds:
            self.master_darks[ccd] = {}
            for exptime in self.obs_run_info['darks'][ccd].keys():
                self.master_darks[ccd][exptime] = {
                    os.path.join(self.obs_run_path, 'darks', 'ccd_1',
                                 exptime, 'DARKcombined_1800.0.fits')}

    def get_master_darks(self, exptime='1800.0'):
        """blah."""
        self.master_darks = {}
        verbose.log_header('Searching master dark files')
        for ccd in self.ccds:
            path_to_master = os.path.join(
                self.obs_run_path, 'darks', ccd, exptime,
                'DARKcombined_{}.fits'.format(exptime))
            if os.path.isfile(path_to_master):
                logging.info('--> [{}] MASTERDARK found at {}\n'.format(
                    ccd, path_to_master))
                self.master_darks[ccd] = path_to_master
            else:
                verbose.missing_master(path_to_master)

    def reduce_lflats(self):
        """blah."""
        verbose.log_header('Reducing long-slit (detector) flats')
        assert self.lflat_idx_file is not None, verbose.missing_idx('lflat_idx')
        self.master_lflats = {}
        for ccd in self.ccds:
            gratings = self.obs_run_info['fflats'][ccd].keys()
            self.master_lflats[ccd] = {}
            for grating in gratings:
                self.master_lflats[ccd][grating] = {}
                for name, fflat_file in self.obs_run_info['fflats'][ccd][grating].items():
                    fflat_exptime = fflat_file['EXPTIME']
                    if fflat_exptime not in self.master_lflats[ccd][grating].keys():
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
                                       '.fits', '_log.txt'),
                                   log=True)
                    self.master_lflats[ccd][grating][fflat_exptime].append(
                        path_to_fflat.replace('.fits', 'red.fits'))
                for exptime in self.master_lflats[ccd][grating].keys():
                    verbose.log_header('Computing MASTERLFLAT [{}] [{}] [{}]\n'.format(
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
                                   output=logmaster,
                                   wdir=os.path.dirname(masterfflat_name),
                                   log=True)
                    self.master_lflats[ccd][grating][exptime] = masterfflat_name
                    logging.info(
                        '---> MASTERLFLAT file saved as %s \n' % masterfflat_name)
                    # Quality control checks
                    QC.check_image(masterfflat_name, save_dir=os.path.join(
                        self.obs_run_path, 'fflats', ccd, grating,
                        'masterlflat_{}.png'.format(exptime)))
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
                if not os.path.isfile(recommended_file):
                    verbose.NoFileError(recommended_file)
                else:
                    with open(recommended_file, 'r') as f:
                        path_to_master = f.readline()
                    if os.path.isfile(path_to_master):
                        logging.info('--> [{}] [{}] MASTERLFLAT found at {}\n'
                                     .format(ccd, grating, path_to_master))
                        self.master_lflats[ccd][grating] = path_to_master
                    else:
                        verbose.missing_master(path_to_master)

    def extract_tramlines(self):
        """blah."""
        verbose.log_header('Starting tramline extraction from fibre flats')
        assert self.fibreflat_idx_file is not None, verbose.missing_idx('fflat_idx')
        self.master_tlm = {}
        for night in self.nights:
            self.master_tlm[night] = {}
            for ccd in self.ccds:
                self.master_tlm[night][ccd] = {}
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
                                           '.fits', '_log.txt'),
                                       log=True)
                        bad_tlm = QC.check_tramline(path_to_fibreflat.replace('.fits', 'tlm.fits'))
                        if not bad_tlm:
                            exptimes.append(float(tram_file['EXPTIME']))
                            names.append(name)
                    recommended_exp_time = kcs.lflat_time[ccd][
                        grating.split('_')[0]]
                    if len(exptimes) > 0:
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
                        self.master_tlm[night][ccd][grating] = best_fibreflat
                        recommended = os.path.join(
                            os.path.dirname(best_fibreflat), 'RECOMMENDED_TRAM')
                        with open(recommended, 'w') as f:
                            f.write(best_fibreflat)
                    else:
                        # TODO
                        pass

    def get_master_tlm(self):
        """blah."""
        self.master_tlm = {}
        verbose.log_header('Searching TRAMLINE MAPS\n')
        for night in self.nights:
            self.master_tlm[night] = {}
            for ccd in self.ccds:
                self.master_tlm[night][ccd] = {}
                gratings = self.obs_run_info[night][ccd].keys()
                for grating in gratings:
                    recommended_file = os.path.join(
                        self.obs_run_path, night, ccd, grating, 'fibreflat',
                        'RECOMMENDED_TRAM')
                    if not os.path.isfile(recommended_file):
                        verbose.NoFileError(recommended_file)
                    else:
                        with open(recommended_file, 'r') as f:
                            path_to_tram = f.readline()
                        if os.path.isfile(path_to_tram):
                            logging.info(
                                '--> [{}] [{}] [{}] Selected Tramline {}\n'.format(
                                    night, ccd, grating, path_to_tram))
                            self.master_tlm[night][ccd][grating] = path_to_tram
                        else:
                            verbose.missing_master(path_to_tram)

    def reduce_arcs(self):
        """blah."""
        verbose.log_header('Reducing calibration arcs')
        assert self.arc_idx_file is not None, verbose.missing_idx('arcs_idx')
        self.master_arcs = {}
        self.check_masters(names=['darks', 'lflats', 'tlm'])
        for night in self.nights:
            self.master_arcs[night] = {}
            for ccd in self.ccds:
                self.master_arcs[night][ccd] = {}
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
                        command_success = aaorun_command(
                            'reduce_arc', path_to_arc,
                            idx_file=self.arc_idx_file,
                            options=['-DARK_FILENAME {}'.format(
                                self.master_darks[ccd]),
                                '-TLMAP_FILENAME {}'.format(
                                    self.master_tlm[night][ccd][grating])],
                            output=path_to_arc.replace('.fits', '_log.txt'),
                            log=True
                        )
                        if command_success == 0:
                            QC.check_image(path_to_arc.replace('.fits', 'red.fits'),
                                           save_dir=os.path.join(
                                               os.path.dirname(path_to_arc),
                                               '{}_arc_{}_{}.png'.format(
                                                   name.split('/')[-1], arc_name, arc_exptime)))
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
                    self.master_arcs[night][ccd][grating] = best_arc
                    recommended = os.path.join(
                        os.path.dirname(best_arc), 'RECOMMENDED_ARC')
                    with open(recommended, 'w') as f:
                        f.write(best_arc)

    def get_master_arcs(self):
        """blah."""
        self.master_arcs = {}
        verbose.log_header('Searching ARCS')
        for night in self.nights:
            self.master_arcs[night] = {}
            for ccd in self.ccds:
                self.master_arcs[night][ccd] = {}
                gratings = self.obs_run_info[night][ccd].keys()
                for grating in gratings:
                    recommended_file = os.path.join(
                        self.obs_run_path, night, ccd, grating, 'arcs',
                        'RECOMMENDED_ARC')
                    if not os.path.isfile(recommended_file):
                        verbose.NoFileError(recommended_file)
                    else:
                        with open(recommended_file, 'r') as f:
                            path_to_arc = f.read().split('\n')[0]
                        if os.path.isfile(path_to_arc):
                            logging.info(
                                '--> [{}] [{}] [{}] ARC found at {}\n'.format(
                                    night, ccd, grating, path_to_arc))
                            self.master_arcs[night][ccd][grating] = path_to_arc
                        else:
                            verbose.missing_master(path_to_arc)

    def reduce_fflats(self):
        """blah."""
        verbose.log_header('Reducing FIBRE FLATS')
        assert self.fibreflat_idx_file is not None, verbose.missing_idx('fflat_idx')
        self.check_masters(names=['darks', 'lflats', 'tlm', 'arcs'])
        self.master_fibreflats = self.master_tlm.copy()
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
                        command_success = aaorun_command(
                            'reduce_fflat', path_to_fflat,
                            idx_file=self.fibreflat_idx_file,
                            options=[
                                '-DARK_FILENAME {}'.format(
                                    self.master_darks[ccd]),
                                '-TLMAP_FILENAME {}'.format(
                                    self.master_tlm[night][ccd][grating]),
                                '-WAVEL_FILENAME {}'.format(
                                    self.master_arcs[night][ccd][grating])],
                            output=path_to_fflat.replace('.fits', '_log.txt'),
                            log=True
                        )
                        if command_success == 0:
                            QC.check_image(path_to_fflat.replace('.fits', 'red.fits'),
                                           save_dir=os.path.join(os.path.dirname(path_to_fflat),
                                                                 '{}.png'.format(name)))
                    self.master_fibreflats[night][ccd][grating] = (
                        self.master_fibreflats[night][ccd][grating].replace(
                            'tlm.fits', 'red.fits'))

    def get_master_fibreflats(self):
        """blah."""
        self.master_fibreflats = {}
        verbose.log_header('Searching fibre flats')
        for night in self.nights:
            self.master_fibreflats[night] = {}
            for ccd in self.ccds:
                self.master_fibreflats[night][ccd] = {}
                gratings = self.obs_run_info[night][ccd].keys()
                for grating in gratings:
                    recommended_file = os.path.join(
                        self.obs_run_path, night, ccd, grating, 'fibreflat',
                        'RECOMMENDED_TRAM')
                    if not os.path.isfile(recommended_file):
                        verbose.NoFileError(recommended_file)
                    else:
                        with open(recommended_file, 'r') as f:
                            path_to_tram = f.readline()
                        path_to_fflat = path_to_tram.replace(
                            'tlm.fits', 'red.fits')
                        if os.path.isfile(path_to_fflat):
                            logging.info(
                                '--> [{}] [{}] [{}] FFLAT found at {}\n'.format(
                                    night, ccd, grating, path_to_fflat))
                            self.master_fibreflats[night][ccd][grating] = path_to_fflat
                        else:
                            verbose.missing_master(path_to_fflat)

    def reduce_object(self):
        """blah."""
        verbose.log_header('Reducing science objects')
        assert self.object_idx_file is not None, verbose.missing_idx('object_idx')
        self.check_masters(names=['darks', 'lflats', 'tlm', 'arcs', 'fflats'])
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
                        command_success = aaorun_command(
                            'reduce_object', path_to_obj,
                            idx_file=self.object_idx_file,
                            options=[
                                '-DARK_FILENAME {}'.format(
                                    self.master_darks[ccd]),
                                '-FFLAT_FILENAME {}'.format(
                                    self.master_fibreflats[night][ccd][
                                        grating]),
                                '-TLMAP_FILENAME {}'.format(
                                    self.master_tlm[night][ccd][grating]),
                                '-WAVEL_FILENAME {}'.format(
                                    self.master_arcs[night][ccd][grating])],
                            output=path_to_obj.replace('.fits', '_log.txt'),
                            log=True
                        )
                        if command_success == 0:
                            QC.check_image(path_to_obj.replace('.fits', 'red.fits'),
                                           save_dir=os.path.join(os.path.dirname(path_to_obj),
                                                                 '{}.png'.format(name)))

    def combine_science_flats(self):
        """Combine all dome/sky flats into a master file for each night."""
        pass


if __name__ == '__main__':
    import time
    tstart = time.time()
    obsrunpath = '/home/pablo/Research/obs_data/HI-KIDS/raw/mar2022_obsruns/obs_run_0'
    redOR = ReduceObsRun(
        obs_run_path=obsrunpath,
        ccds=['ccd_1', 'ccd_2'],
        dark_idx='koala_dark.idx',
        fibreflat_idx='koala_dark.idx',
        lflat_idx='koala_dark.idx',
        arcs_idx='koala_arcs.idx',
        object_idx='koala_reduce.idx')

    # redOR.reduce_darks()
    redOR.reduce_lflats()

    #

    # redOR.get_master_darks()
    # redOR.get_master_lflats()
    # redOR.get_master_tlm()
    # redOR.get_arcs()
    # redOR.get_fibreflats()
    # redOR.extract_tramlines()

    # redOR.reduce_arcs()
    # redOR.reduce_fflats()
    #redOR.reduce_object()

    tend = time.time()
    print('\n\n ### Elapsed time (hrs): ', (tend - tstart) / 3600)

# Mr Krtxo \(ﾟ▽ﾟ)/
