#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains...

Log of (some) changes
08-08-22:

"""

import yaml
import os
from .py2dfdr import aaorun_command
from . import reduction_QC as QC
from . import koala_cheatsheet as kcs
import numpy as np
import logging
import datetime
from . import verbose


class ReduceObsRun(object):
    """Reduce KOALA observing runs.

    Description
    -----------
    This class provides performs all the necessary reduction steps and quality control checks to reduce an observing run
    from KOALA data.

    Attributes
    ----------
    - obs_run_info: (dict) Description of the contents within the OR. When providing the path to a given OR, the
    obs_run_info.yml must be on the root path.
    - obs_run_path: (str) Path to the root of the OR.
    - ccds: (list) List of CCDs detectors to reduce (e.g. "[ccd_1, ccd_2]")
    - dark_idx_file: (str) Name of the .idx file used to reduce dark files.
    - lflat_idx_file: (str) Name of the .idx file used to reduce long-slit flat files.
    - fibreflat_idx_file: (str) Name of the .idx file used to reduce fibre flat files (including tramline mapping).
    - dark_idx_file: (str) Name of the .idx file used to reduce dark files.
    - nights: (list) List containing all the observing nights within the OR.
    - master_bias: (dict) Dictionary containing the information of master bias.
    - master_darks: (dict) Dictionary containing the information of master darks (one for each exp. time).
    - master_lflats: (dict) Dictionary containing the information of master long-slit flats (
    one for each exp. time and grating).
    - master_tlm: (dict) Dictionary containing the information of tramline maps (one for each night, ccd and grating)
    - master_arcs: (dict) Dictionary containing the information of arc lamps file (one for each night, ccd and grating)
    - master_fibreflats: (dict) Dictionary containing the information of fibreflats (one for each night, ccd and
    grating)
    - sat_fraction: (float, default=0.3) Maximum fraction of saturated pixels for preforming the reduction. Files with
    values above this limit will be flagged as "SATURATED" and no reduction will be performed.
    - sat_level: (float, default=65500.0) Number of counts for saturated pixels.
    - reject_names: (str list, default=['FOCUS']) List containing names/keywords of files to reject during data
    reduction based on the fits header.
    - verb: (bool, default=True) Whether to print reduction steps or only saving them on the log file.

    Methods
    -------
    - load_obs_run
    - check_masters
    - reject_saturated
    - reject_from_name
    - reduce_bias
    - get_master_bias
    - reduce_darks
    - get_master_darks
    - reduce_lflats
    - get_master_lflats
    - extract_tramlines
    - get_master_tlm
    - reduce_arcs
    - get_master_arcs
    - reduce_fflats
    - get_master_fflats
    - reduce_object
    - combine_science_data

    Example
    -------
    tstart = time.time()
    # Path to the observing run
    obsrunpath = '/home/pablo/Research/obs_data/HI-KIDS/raw/mar2022_obsruns/obs_run_0'
    # Instantiate the ReduceObsRun object
    redOR = red_obs.ReduceObsRun(obs_run_path=obsrunpath, ccds=['ccd_1', 'ccd_2'],
                                 dark_idx='koala_dark.idx',
                                 lflat_idx='koala_dark.idx',
                                 fibreflat_idx='koala_fflat.idx',
                                 arcs_idx='koala_arcs.idx',
                                 object_idx='koala_reduce.idx',
                                 # Values above or equal to 65500 will be considered saturated pixels
                                 sat_level=65500.0,
                                 # When the fraction of saturated pixels is >= sat_fraction the file is not reduced
                                 sat_fraction=0.5)
    # Start the reduction sequence
    redOR.reduce_darks(timeout=300)  # If the aaorun command takes more than "timeout", the process is skipped.
    # redOR.get_master_darks()
    redOR.reduce_lflats(timeout=300)
    # redOR.get_master_lflats()
    redOR.extract_tramlines(timeout=300)
    # redOR.get_master_tlm()
    # redOR.reduce_arcs(timeout=300)
    # redOR.get_arcs()
    # redOR.reduce_fflats(timeout=300)
    # redOR.get_fibreflats()
    redOR.reduce_object(timeout=900)
    tend = time.time()
    print('\n\n ### Elapsed time (hrs): ', (tend - tstart) / 3600)
    # Good luck ;)
    # Mr Krtxo \(ﾟ▽ﾟ)/
    """

    def __init__(self, obs_run_path, verb=True, **kwargs):
        """
        Observing Run Reduction constructor

        Input params
        ------------
        - obs_run_path: (str) Path to OR directory
        - verb: (bool, optional, default=True) If True, all the process output will be printed as recorded on the Log
        file.
        - kwargs:
            - sat_frac
            - sat_level
            - reject_names
            - dark_idx
            - lflat_idx
            - fibreflat_idx
            - arcs_idx
            - object_idx
        """
        self.obs_run_info = None
        self.nights = None
        self.master_bias = None
        self.master_darks = None
        self.master_lflats = None
        self.master_tlm = None
        self.master_arcs = None
        self.master_fibreflats = None

        # Data rejection
        self.sat_fraction = kwargs.get('sat_frac', 0.3)
        self.sat_level = kwargs.get('sat_level', 65500.)
        self.reject_names = kwargs.get('reject_names', ['FOCUS'])

        self.obs_run_path = obs_run_path
        # Initialise logging file
        logging.basicConfig(
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%m/%d/%Y %I:%M:%S %p',
            filename=os.path.join(self.obs_run_path, 'OR_reduction.log'),
            # handlers=[
            #    logging.FileHandler("OR_reduction.log"),
            #    logging.StreamHandler()],
            level=logging.INFO)
        if verb:
            logging.getLogger().addHandler(logging.StreamHandler())
        verbose.log_header('[OBSRUN] Initialising OR reduction process at:\n   {} \n'.format(
            self.obs_run_path))
        logging.info(datetime.datetime.now().strftime("%c"))
        # CCD detector to reduce
        self.ccds = kwargs.get('ccds', None)
        if self.ccds is None:
            logging.error('[OBSRUN] ERROR: CCDs not provided [ccd_1, ccd_2] \n')
            raise NameError('[OBSRUN] ERROR: CCDs not provided [ccd_1, ccd_2]')
        else:
            logging.info('CCDs data to reduce: {}'.format(', '.join(self.ccds)))
        # Parameter files for 2dfdr
        logging.info('[OBSRUN] Setting configuration files for 2dfdr')
        self.dark_idx_file = kwargs.get('dark_idx', None)
        logging.info('[OBSRUN] DARK configuration file: %s' % self.dark_idx_file)
        self.lflat_idx_file = kwargs.get('lflat_idx', None)
        logging.info('[OBSRUN] LFLAT configuration file: %s' % self.lflat_idx_file)
        self.fibreflat_idx_file = kwargs.get('fibreflat_idx', None)
        logging.info('[OBSRUN] FFLAT configuration file: %s' % self.fibreflat_idx_file)
        self.arc_idx_file = kwargs.get('arcs_idx', None)
        logging.info('[OBSRUN] ARCS configuration file: %s' % self.arc_idx_file)
        self.object_idx_file = kwargs.get('object_idx', None)
        logging.info('[OBSRUN] OBJECT configuration file: %s' % self.object_idx_file)

        # Load yml file containing data description
        self.load_obs_run_yml(**kwargs)

    def load_obs_run_yml(self, **kwargs):
        """Load the Observing Run yaml file."""
        logging.info('[OBSRUN] · Loading yml file containing the OR data\n')
        with open(os.path.join(self.obs_run_path, "obs_run_info.yml"),
                  "r") as stream:
            try:
                self.obs_run_info = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                logging.error('[OBSRUN] · ERROR: Unable to load yml file\n')
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
                if masters[name] is None:
                    get_method[name]()

    # REJECTION METHODS ------------------------------------------------------------------------------------------------
    def reject_saturated_image(self, path):
        """Reject object based on fraction of saturated pixels."""
        sat_frac = QC.check_saturated(path, sat_level=self.sat_level, log=True)
        if sat_frac >= self.sat_fraction:
            return True
        else:
            return False

    def reject_from_name(self, name):
        """Reject object based on header name."""
        for rej_name in self.reject_names:
            if name.find(rej_name) >= 0:
                return True
        return False

    @staticmethod
    def reject_binned(path):
        """Reject data on binning mode."""
        val = QC.get_keyword(path, keyword='WINDOW', hdu_index=0)
        if val is not None:
            if val.find('BIN') >= 0:
                return True
            else:
                return False
        else:
            return False

    # REDUCTION METHODS ------------------------------------------------------------------------------------------------
    def reduce_bias(self):
        """blah."""
        verbose.log_header('bias')
        for ccd in self.ccds:
            for name, bias_file in self.obs_run_info['bias'][ccd].items():
                logging.info('BIAS: %s\n' % bias_file['PATH'])

    def get_master_bias(self):
        pass

    def reduce_darks(self, timeout=900):
        """blah."""
        verbose.log_header('Reducing darks')
        if self.dark_idx_file is None:
            verbose.missing_idx('dark_idx')
        for ccd in self.ccds:
            for exptime in self.obs_run_info['darks'][ccd].keys():
                all_darks = []
                for name, dark_file in self.obs_run_info['darks'][ccd][exptime].items():
                    logging.info('\n[OBSRUN] ·[{}] [{}] DARK: {}\n'.format(
                        ccd, exptime,
                        dark_file['PATH']))

                    path_to_dark = os.path.join(
                        self.obs_run_path, 'darks', ccd, exptime,
                        dark_file['PATH'])
                    # CALL AAORUN
                    aaorun_command('reduce_dark', path_to_dark,
                                   idx_file=self.dark_idx_file,
                                   output=path_to_dark.replace(
                                       '.fits', '_log.txt'),
                                   timeout=timeout)
                    all_darks.append(path_to_dark.replace('.fits', 'red.fits'))
                logging.info(
                    '[OBSRUN] · Computing MASTERDARK [{}] [{}]\n'.format(ccd, exptime))
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
                    '[OBSRUN] MASTERDARK file saved as %s \n' % masterdark_name)
                logging.info('[OBSRUN] MASTERDARK log saved as %s \n' % logmaster)
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
                                 exptime, 'DARKcombined_{}.fits'.format(exptime))}

    def get_master_darks(self, exptime=None):
        """blah."""
        self.master_darks = {}
        verbose.log_header('[OBSRUN] · Searching master dark files')
        for ccd in self.ccds:
            if exptime is None:
                best_dark = np.argmax(np.array(list(self.obs_run_info['darks'][ccd].keys()), dtype=np.float))
                exptime = list(self.obs_run_info['darks'][ccd].keys())[best_dark]
            path_to_master = os.path.join(
                self.obs_run_path, 'darks', ccd, exptime,
                'DARKcombined_{}.fits'.format(exptime))
            if os.path.isfile(path_to_master):
                logging.info('[OBSRUN] [{}] MASTERDARK found at {}'.format(
                    ccd, path_to_master))
                self.master_darks[ccd] = path_to_master
            else:
                verbose.missing_master(path_to_master)

    def reduce_lflats(self, timeout=900):
        """blah."""
        verbose.log_header('Reducing long-slit (detector) flats')
        if self.lflat_idx_file is None:
            verbose.missing_idx('lflat_idx')
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
                    logging.info('\n[OBSRUN] ·[{}] [{}] LFLAT: {}'.format(
                        ccd, grating, fflat_file['PATH']))
                    path_to_fflat = os.path.join(
                        self.obs_run_path, 'fflats', ccd, grating,
                        fflat_file['PATH'])
                    # CALL TO AAORUN
                    aaorun_command('reduce_lflat', path_to_fflat,
                                   idx_file=self.lflat_idx_file,
                                   output=path_to_fflat.replace(
                                       '.fits', '_log.txt'),
                                   log=True,
                                   timeout=timeout)
                    self.master_lflats[ccd][grating][fflat_exptime].append(
                        path_to_fflat.replace('.fits', 'red.fits'))

                exposure_times = list(self.master_lflats[ccd][grating].keys())
                for exptime in exposure_times:
                    verbose.log_header('[OBSRUN] Computing MASTERLFLAT [{}] [{}] [{}]\n'.format(
                        ccd, grating, exptime))
                    if len(self.master_lflats[ccd][grating][exptime]) > 1:
                        fflats_to_combine = ' '.join(
                            self.master_lflats[ccd][grating][exptime])
                    else:
                        del self.master_lflats[ccd][grating][exptime]
                        logging.warning(
                            '[OBSRUN] Only one file available! Skipping combination')
                        continue
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
                        '[OBSRUN] MASTERLFLAT file saved as\n  %s' % masterfflat_name)
                    # Quality control checks
                    QC.check_image(masterfflat_name, save_dir=os.path.join(
                        self.obs_run_path, 'fflats', ccd, grating,
                        'masterlflat_{}.png'.format(exptime)))
                # Best exposure time
                logging.info('\n[OBSRUN] ·[{}] [{}] LFLAT: Selecting best master LFLAT'.format(
                    ccd, grating))
                recommended_exp_time = kcs.lflat_time[ccd][grating.split('_')[0]]
                time_keys = list(self.master_lflats[ccd][grating].keys())
                times = np.array(time_keys, dtype=float)
                best = np.argmin(np.abs(times - recommended_exp_time))
                logging.info('\n[OBSRUN] ·[{}] [{}] LFLAT: Recommended exp. time {} \n Best available: {}, {}'.format(
                    ccd, grating, recommended_exp_time, times[best], self.master_lflats[ccd][grating][time_keys[best]]))
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
        logging.info('[OBSRUN] · Searching master dark files')
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
                        logging.info('[OBSRUN] [{}] [{}] MASTERLFLAT found at\n  {}'
                                     .format(ccd, grating, path_to_master))
                        self.master_lflats[ccd][grating] = path_to_master
                    else:
                        verbose.missing_master(path_to_master)

    def extract_tramlines(self, timeout=900):
        """blah."""
        verbose.log_header('Starting tramline extraction from fibre flats')
        if self.fibreflat_idx_file is None:
            verbose.missing_idx('fflat_idx')
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
                        logging.info('\n[OBSRUN] · [{}] [{}] [{}]'.format(night, ccd,
                                                                          grating)
                                     + ' TRAM: %s' % tram_file['PATH'])
                        path_to_fibreflat = os.path.join(
                            self.obs_run_path, night, ccd, grating,
                            'fibreflat', tram_file['PATH'])
                        if self.reject_saturated_image(path_to_fibreflat):
                            logging.info('[OBSRUN] · [{}] [{}] [{}]'.format(night, ccd,
                                                                            grating)
                                         + ' TRAM: SATURATED! SKIPPING REDUCTION.')
                            continue
                        # CALL TO AAORUN
                        command_success = aaorun_command('make_tlm', path_to_fibreflat,
                                                         idx_file=self.fibreflat_idx_file,
                                                         output=path_to_fibreflat.replace(
                                                         '.fits', '_log.txt'),
                                                         log=True,
                                                         timeout=timeout)
                        if command_success == 0:
                            bad_tlm = QC.check_tramline(path_to_fibreflat.replace('.fits', 'tlm.fits'))
                        else:
                            bad_tlm = True
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
                        logging.info('[OBSRUN] · [{}] [{}] [{}]'.format(night, ccd, grating)
                                     + ' BEST RECOMMENDED TRAM:\n  %s'
                                     % best_fibreflat)
                        self.master_tlm[night][ccd][grating] = best_fibreflat
                        recommended = os.path.join(
                            os.path.dirname(best_fibreflat), 'RECOMMENDED_TRAM')
                        with open(recommended, 'w') as f:
                            f.write(best_fibreflat)
                    else:
                        logging.warning('[OBSRUN] · [{}] [{}] [{}]'.format(night, ccd, grating)
                                        + 'WARNING: NO TRAMLINE AVAILABLE.')

    def get_master_tlm(self):
        """blah."""
        self.master_tlm = {}
        verbose.log_header('[OBSRUN] · Searching TRAMLINE MAPS')
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
                                '[OBSRUN] [{}] [{}] [{}] Selected Tramline\n  {}'.format(
                                    night, ccd, grating, path_to_tram))
                            self.master_tlm[night][ccd][grating] = path_to_tram
                        else:
                            verbose.missing_master(path_to_tram)

    def reduce_arcs(self, timeout=900):
        """blah."""
        verbose.log_header('[OBSRUN] · Reducing calibration arcs')
        if self.arc_idx_file is None:
            verbose.missing_idx('arcs_idx')
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
                        logging.info('\n[OBSRUN] · [{}] [{}] [{}] [{}] [{}]'.format(
                            night, ccd, grating, arc_name, arc_exptime)
                                     + ' ARC: %s' % arc_file['PATH'])
                        path_to_arc = os.path.join(
                            self.obs_run_path, night, ccd, grating, 'arcs',
                            arc_file['PATH'])
                        if self.reject_saturated_image(path_to_arc):
                            logging.info('[OBSRUN] · [{}] [{}] [{}]'.format(night, ccd,
                                                                            grating)
                                         + ' ARC SATURATED! SKIPPING REDUCTION.')
                            continue
                        # CALL TO AAORUN
                        command_success = aaorun_command(
                            'reduce_arc', path_to_arc,
                            idx_file=self.arc_idx_file,
                            options=['-DARK_FILENAME {}'.format(
                                self.master_darks[ccd]),
                                '-TLMAP_FILENAME {}'.format(
                                    self.master_tlm[night][ccd][grating])],
                            output=path_to_arc.replace('.fits', '_log.txt'),
                            log=True,
                            timeout=timeout
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
                    if selected_arcs.any():
                        best = np.nanargmin(
                            np.abs(exptimes[selected_arcs] - rec_exptime))
                        best_arc = os.path.join(
                            self.obs_run_path, night, ccd,
                            grating, 'arcs',
                            self.obs_run_info[night][ccd][grating]['arcs'][
                                names[best]]['PATH'].replace('.fits', 'red.fits'))
                        logging.info(
                            '[OBSRUN] · [{}] [{}] [{}]'.format(night, ccd, grating)
                            + ' BEST RECOMMENDED ARC:\n  %s' % best_arc)
                        self.master_arcs[night][ccd][grating] = best_arc
                        recommended = os.path.join(
                            os.path.dirname(best_arc), 'RECOMMENDED_ARC')
                        with open(recommended, 'w') as f:
                            f.write(best_arc)
                    else:
                        logging.warning(
                            '[OBSRUN] · [{}] [{}] [{}]'.format(night, ccd, grating)
                            + ' WARNING: NO ARC SELECTED!')

    def get_master_arcs(self):
        """blah."""
        self.master_arcs = {}
        verbose.log_header('[OBSRUN] · Searching ARCS')
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
                                '[OBSRUN] [{}] [{}] [{}] ARC found at\n  {}'.format(
                                    night, ccd, grating, path_to_arc))
                            self.master_arcs[night][ccd][grating] = path_to_arc
                        else:
                            verbose.missing_master(path_to_arc)

    def reduce_fflats(self, timeout=900):
        """blah."""
        verbose.log_header('[OBSRUN] · Reducing FIBRE FLATS')
        if self.fibreflat_idx_file is None:
            verbose.missing_idx('fflat_idx')
        self.check_masters(names=['darks', 'lflats', 'tlm', 'arcs'])
        self.master_fibreflats = self.master_tlm.copy()
        # TODO: This should only reduce the master fflat.
        for night in self.nights:
            for ccd in self.ccds:
                gratings = self.obs_run_info[night][ccd].keys()
                for grating in gratings:
                    for name, object_file in self.obs_run_info[night][ccd][
                            grating]['fibreflat'].items():
                        exptime = object_file['EXPTIME']
                        logging.info('\n[OBSRUN] · [{}] [{}] [{}] [{}] '.format(
                            night, ccd, grating, exptime)
                                     + ' FFLAT: %s' % object_file['PATH'])
                        path_to_fflat = os.path.join(self.obs_run_path, night,
                                                     ccd, grating, 'fibreflat',
                                                     object_file['PATH'])
                        if self.reject_saturated_image(path_to_fflat):
                            logging.info('[OBSRUN] · [{}] [{}] [{}]'.format(night, ccd,
                                                                            grating)
                                         + ' FFLAT SATURATED! SKIPPING REDUCTION.')
                            continue
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
                            log=True,
                            timeout=timeout
                        )
                        if command_success == 0:
                            QC.check_image(path_to_fflat.replace('.fits', 'red.fits'),
                                           save_dir=os.path.join(os.path.dirname(path_to_fflat),
                                                                 '{}.png'.format(name)))
                        else:
                            logging.warning('[OBSRUN] · [{}] [{}] [{}]'.format(night, ccd,
                                                                               grating)
                                            + ' Unsuccessful operation')
                    self.master_fibreflats[night][ccd][grating] = (
                        self.master_fibreflats[night][ccd][grating].replace(
                            'tlm.fits', 'red.fits'))

    def get_master_fibreflats(self):
        """blah."""
        self.master_fibreflats = {}
        verbose.log_header('[OBSRUN] · Searching fibre flats')
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
                                '[OBSRUN] [{}] [{}] [{}] FFLAT found at:\n  {}'.format(
                                    night, ccd, grating, path_to_fflat))
                            self.master_fibreflats[night][ccd][grating] = path_to_fflat
                        else:
                            verbose.missing_master(path_to_fflat)

    def reduce_object(self, timeout=900):
        """blah."""
        verbose.log_header('[OBSRUN] · Reducing science objects')
        if self.object_idx_file is None:
            verbose.missing_idx('object_idx')
        self.check_masters(names=['darks', 'lflats', 'tlm', 'arcs', 'fflats'])
        for night in self.nights:
            for ccd in self.ccds:
                gratings = self.obs_run_info[night][ccd].keys()
                for grating in gratings:
                    # Dictionary containing reduction tags
                    object_flags = {}
                    for name, object_file in self.obs_run_info[night][ccd][
                            grating]['sci'].items():
                        obj_name = object_file['NAME']
                        exptime = object_file['EXPTIME']
                        object_flags[name] = {}
                        object_flags[name]['NAME'] = obj_name
                        object_flags[name]['PATH'] = object_file['PATH']
                        logging.info('\n[OBSRUN] · [{}] [{}] [{}] [{}] [{}]'.format(
                            night, ccd, grating, obj_name, exptime)
                                     + ' OBJECT: %s' % object_file['PATH'])
                        path_to_obj = os.path.join(
                            self.obs_run_path, night, ccd, grating, 'sci',
                            object_file['PATH'])
                        # Data rejection
                        if self.reject_from_name(obj_name):
                            object_flags[name]['FLAG'] = 'NAMEREJ'
                            logging.info('[OBSRUN] · [{}] [{}] [{}]'.format(night, ccd,
                                                                            grating)
                                         + ' OBJECT NAME-REJECTED! SKIPPING REDUCTION.')
                            continue
                        if self.reject_binned(path_to_obj):
                            object_flags[name]['FLAG'] = 'BINNING_CONFIG'
                            logging.info('[OBSRUN] · [{}] [{}] [{}]'.format(night, ccd,
                                                                            grating)
                                         + ' OBJECT IN BINNING CONFIGURATION. SKIPPING REDUCTION.')
                            continue
                        if self.reject_saturated_image(path_to_obj):
                            object_flags[name]['FLAG'] = 'SATURATED'
                            logging.info('[OBSRUN] · [{}] [{}] [{}]'.format(night, ccd,
                                                                            grating)
                                         + ' OBJECT SATURATED! SKIPPING REDUCTION.')
                            continue
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
                            log=True,
                            timeout=timeout
                        )
                        if command_success == 0:
                            object_flags[name]['FLAG'] = 'OK'
                            QC.check_image(path_to_obj.replace('.fits', 'red.fits'),
                                           save_dir=os.path.join(os.path.dirname(path_to_obj),
                                                                 '{}.png'.format(name)),
                                           title=obj_name)
                        else:
                            object_flags[name]['FLAG'] = 'AAORUNFAIL'
                            logging.warning('[OBSRUN] · [{}] [{}] [{}]'.format(night, ccd,
                                                                               grating)
                                            + 'WARNING: Unsuccessful reduction')
                    with open(os.path.join(self.obs_run_path, night, ccd, grating, 'sci',
                                           'REDUCTION_FLAGS.yml'), 'w') as outfile:
                        yaml.dump(object_flags, outfile, default_flow_style=False)

    # Extra functions
    def combine_science_data(self, keyword=None):
        """Combine all files within each night sharing a common name."""
        pass


# Mr Krtxo \(ﾟ▽ﾟ)/
