#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import subprocess
import sys
import numpy as np
from astropy.io import fits
from glob import glob
from astropy.time import Time
import shutil
import yaml
import logging


class KOALA_archive(object):
    """This class allows the user to classify and extract all the information from a set of observing runs.

    Attributes
    ----------
    - nights: (dict) Collection ob observing nights with file information.
    - dates: (array) Observing nights julian dates.
    - observing_runs: (dict) Collection of observing runs with file information.
    - root_path: (str) Path to the directory containing all the individual observing nights.

    Methods
    -------
    - find_nights
    - get_night_info
    - get_content
    - filter_night_files
    - create_observing_runs
    - create_obs_run_info
    - make_dir

    Example
    -------
    # Create an archive classifier
    archive = KOALA_archive(path_to_data)
    # Find all the nigths within 'path_to_data'
    archive.find_nights()
    # Extract the information of each file for every night
    archive.get_night_info()
    # Create a directory containing each observing run using symbolic links
    archive.create_reduction_folder(output_path=new_working_dir_path, hard=False)
    """

    def __init__(self, root_path):
        self.nights = None
        self.dates = None
        self.observing_runs = None
        # ------------------------
        self.root_path = root_path
        if self.root_path[-1] != '/':
            self.root_path += '/'

    def find_nights(self):
        """blah."""
        nights_path = np.sort(glob(self.root_path + '*'))
        # Remove any file that could be within the root folder
        nights_path = list(filter(os.path.isdir, nights_path))
        # Check if there is any night
        if len(nights_path) == 0:
            raise NameError(
                'Root path:\n{}\n contains no files'.format(self.root_path))
        else:
            print('[Archiving] {} nights found'.format(len(nights_path)))
            nights = {}
            self.dates = np.zeros(len(nights_path))
            for i, night_dir in enumerate(nights_path):
                folder = night_dir.split('/')[-1]
                julian_date = Time(
                    folder[:4] + '-' + folder[4:6] + '-'+folder[6:]).jd
                self.dates[i] = julian_date
                nights[folder] = {'JD': julian_date}
            sortpos = np.argsort(self.dates)

            self.nights = nights

    def get_night_info(self):
        """blah..."""
        for night in self.nights.keys():
            night_path = os.path.join(self.root_path, night)
            ccd_folders = glob(os.path.join(night_path, '*'))
            ccd_folders = list(filter(os.path.isdir, ccd_folders))
            self.nights[night]['content'] = {}
            for ccd in ccd_folders:
                content = self.get_content(ccd)
                self.nights[night]['content'][ccd.split('/')[-1]] = content

    @staticmethod
    def get_content(path):
        """
        Return a dictionary containing the description of each file within a given path folder.

        params:
        ------
        - path (str): path to the data

        returns:
        --------
         - files: (dict)
        """
        # change the working directory
        os.chdir(path)
        print('Searching files at: ', os.getcwd())
        # Call aaorun list for selecting the elements within the folder
        aaorun_list = os.system('aaorun list > content.txt')
        if aaorun_list != 0:
            raise NameError('ERROR: aaorun list command did not work!')
        # Read the ouput from aaorun and create the "files" dictionary
        with open('content.txt', 'r') as f:
            content = f.read()
            content = content.split('\n')
            print('Header line --> ', content[0])
            files_prefix = content[0][
                content[0].find(': ')+2: content[0].find(' (')]
            print(' > Files prefix *{}*'.format(files_prefix))
            content = content[1:-1]
            files = {}
            for file in range(len(content)):
                data = content[file].split(' ')
                file_number = data[0].split(files_prefix)[1].split('.fits')[0]
                files[file_number] = {'PATH': data[0], 'DESC': data[1]}
                if files[file_number]['DESC'] == 'BIAS':
                    continue
                with fits.open(data[0]) as hdul:

                    if files[file_number]['DESC'] == 'DARK':
                        exptime = hdul[0].header['EXPOSED']
                        files[file_number]['EXPTIME'] = str(exptime)
                        continue
                    else:
                        files[file_number]['GRATID'] = hdul[0].header['GRATID']
                        files[file_number]['GRATWAVE'] = hdul[0].header['LAMBDBR']

                    if data[1] == 'MFOBJECT':
                        name = hdul[0].header['OBJECT'].upper()
                        exptime = hdul[0].header['EXPOSED']
                        files[file_number]['NAME'] = name
                        files[file_number]['EXPTIME'] = str(exptime)
                    elif data[1] == 'MFARC':
                        # Save arc lamp name (e.g. (CuAr)) removing the parenthesis
                        files[file_number]['ARCNAME'] = data[2][1:-1]
                        exptime = hdul[0].header['EXPOSED']
                        files[file_number]['EXPTIME'] =  str(exptime)
                    elif data[1] == 'LFLAT':
                        exptime = hdul[0].header['EXPOSED']
                        files[file_number]['EXPTIME'] =  str(exptime)
                    elif data[1] == 'MFFFF':
                        exptime = hdul[0].header['EXPOSED']
                        files[file_number]['EXPTIME'] =  str(exptime)
                hdul.close()
        f.close()
        print(' > Found {} files'.format(len(files)))
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        return files

    def filter_night_files(self):
        pass

    def create_observing_runs(self):
        """blah."""
        print(' [Archiving] Grouping nights into observing runs')
        if len(self.dates) > 1:
            diffdates = np.diff(self.dates)
            binnumber = np.arange(0, len(self.dates), 1)
            for i in range(diffdates.size):
                if diffdates[i] == 1:
                    binnumber[i+1] = binnumber[i]
            self.observing_run_id = np.unique(binnumber)
            self.night_run_id = binnumber
            self.nights_keys = np.array(list(self.nights.keys()))
            for i, night in enumerate(self.nights.keys()):
                self.nights[night]['obsrun'] = binnumber[i]
        else:
            binnumber = [0]
            self.observing_run_id = binnumber
            self.night_run_id = binnumber
            self.nights_keys = np.array(list(self.nights.keys()))
            for i, night in enumerate(self.nights_keys):
                self.nights[night]['obsrun'] = binnumber[i]

    def create_reduction_folder(self, output_path=None, hard=True):
        """blah..."""
        # Create root dir for reduction
        if output_path is None:
            output_path = self.root_path
        if hard:
            copier = shutil.copy
        else:
            copier = os.symlink

        reduce_path = os.path.join(output_path)
        self.make_dir(reduce_path)
        # Create observing runs folders
        self.observing_runs = {}
        for obs_run in self.observing_run_id:
            self.observing_runs[obs_run] = {}
            obs_run_path = os.path.join(reduce_path,
                                        'obs_run_{}'.format(obs_run))
            self.make_dir(obs_run_path)
            # Create common files dir
            self.make_dir(os.path.join(obs_run_path, 'darks'))
            self.make_dir(os.path.join(obs_run_path, 'darks', 'ccd_1'))
            self.make_dir(os.path.join(obs_run_path, 'darks', 'ccd_2'))
            self.observing_runs[obs_run]['darks'] = {'ccd_1': {}, 'ccd_2': {}}

            self.make_dir(os.path.join(obs_run_path, 'bias'))
            self.make_dir(os.path.join(obs_run_path, 'bias', 'ccd_1'))
            self.make_dir(os.path.join(obs_run_path, 'bias', 'ccd_2'))
            self.observing_runs[obs_run]['bias'] = {'ccd_1': {}, 'ccd_2': {}}

            self.make_dir(os.path.join(obs_run_path, 'fflats'))
            self.make_dir(os.path.join(obs_run_path, 'fflats', 'ccd_1'))
            self.make_dir(os.path.join(obs_run_path, 'fflats', 'ccd_2'))
            self.observing_runs[obs_run]['fflats'] = {'ccd_1': {}, 'ccd_2': {}}

            obs_run_nights_id = (
                self.nights_keys[np.where(self.night_run_id == obs_run)[0]])
            for nightid in obs_run_nights_id:
                # Night data
                night = self.nights[nightid]
                # Path to night data
                path_to_night = os.path.join(self.root_path, nightid)
                # New path to obs_run night data
                obs_run_night_path = os.path.join(
                    obs_run_path, 'night_' + nightid)
                # Night-specific data
                self.make_dir(obs_run_night_path)
                self.make_dir(os.path.join(obs_run_night_path, 'ccd_1'))
                self.make_dir(os.path.join(obs_run_night_path, 'ccd_2'))
                self.observing_runs[obs_run]['night_' + nightid] = {
                    'ccd_1': {}, 'ccd_2': {}}

                for ccd in night['content'].keys():
                    path_to_night_ccd = os.path.join(path_to_night, ccd)
                    # Loop over all the files in the night folder
                    for fid, file in night['content'][ccd].items():
                        if file['DESC'] == 'DARK':
                            exptime = file['EXPTIME']
                            if not os.path.isdir(os.path.join(
                                    obs_run_path, 'darks', ccd, exptime)):
                                self.make_dir(os.path.join(
                                        obs_run_path, 'darks', ccd, exptime))
                                self.observing_runs[
                                    obs_run]['darks'][ccd][exptime] = {}

                            copier(os.path.join(path_to_night_ccd,
                                                file['PATH']),
                                   os.path.join(obs_run_path, 'darks',
                                                ccd, exptime, file['PATH']))
                            self.observing_runs[
                                obs_run]['darks'][ccd][exptime]['_'.join((nightid, fid))] = (file)
                        elif file['DESC'] == 'BIAS':
                            copier(os.path.join(path_to_night_ccd,
                                                file['PATH']),
                                   os.path.join(obs_run_path, 'bias',
                                                ccd, file['PATH']))
                            self.observing_runs[obs_run]['bias'][ccd]['_'.join((nightid, fid))] = (
                                file)
                        elif file['DESC'] == 'LFLAT':
                            grating = file['GRATID'] + '_{:.0f}'.format(
                                file['GRATWAVE'])
                            if not os.path.isdir(os.path.join(
                                    obs_run_path, 'fflats', ccd, grating)):
                                self.make_dir(os.path.join(
                                        obs_run_path, 'fflats', ccd,
                                        grating))
                                self.observing_runs[obs_run]['fflats'][ccd][
                                    grating] = {}
                            copier(os.path.join(path_to_night_ccd,
                                                file['PATH']),
                                   os.path.join(obs_run_path, 'fflats',
                                                ccd, grating,
                                                file['PATH']))
                            self.observing_runs[obs_run]['fflats'][ccd][
                                grating]['_'.join((nightid, fid))] = (file)
                        else:
                            grating = file['GRATID'] + '_{:.0f}'.format(
                                file['GRATWAVE'])

                            if not os.path.isdir(os.path.join(
                                    obs_run_night_path, ccd, grating)):
                                self.make_dir(os.path.join(
                                    obs_run_night_path, ccd, grating))
                                self.make_dir(os.path.join(
                                    obs_run_night_path, ccd, grating,
                                    'sci'))
                                self.make_dir(os.path.join(
                                    obs_run_night_path, ccd, grating,
                                    'fibreflat'))
                                self.make_dir(os.path.join(
                                    obs_run_night_path, ccd, grating,
                                    'arcs'))

                                self.observing_runs[obs_run]['night_' + nightid][ccd][
                                    grating] = {'sci': {}, 'fibreflat': {},
                                                'arcs': {}}

                            if file['DESC'] == 'MFFFF':
                                copier(os.path.join(path_to_night_ccd,
                                                    file['PATH']),
                                       os.path.join(obs_run_night_path,
                                                    ccd, grating, 'fibreflat',
                                                    file['PATH']))
                                self.observing_runs[obs_run][
                                    'night_' + nightid][ccd][
                                    grating]['fibreflat']['_'.join((nightid, fid))] = file
                            elif file['DESC'] == 'MFOBJECT':
                                copier(os.path.join(path_to_night_ccd,
                                                    file['PATH']),
                                       os.path.join(obs_run_night_path,
                                                    ccd, grating, 'sci',
                                                    file['PATH']))
                                self.observing_runs[obs_run][
                                    'night_' + nightid][ccd][
                                    grating]['sci']['_'.join((nightid, fid))] = file
                            elif file['DESC'] == 'MFARC':
                                copier(os.path.join(path_to_night_ccd,
                                                    file['PATH']),
                                       os.path.join(obs_run_night_path,
                                                    ccd, grating, 'arcs',
                                                    file['PATH']))
                                self.observing_runs[obs_run][
                                    'night_' + nightid][ccd][
                                    grating]['arcs']['_'.join((nightid, fid))] = file
            self.create_obs_run_info(self.observing_runs[obs_run],
                                     savepath=os.path.join(obs_run_path,
                                                           'obs_run_info'))

    def create_obs_run_info(self, obsrun, savepath):
        """Create a yml file to store the obs run data."""
        with open(savepath + '.yml', 'w') as outfile:
            yaml.dump(obsrun, outfile, default_flow_style=False)

    def make_dir(self, path):
        """Create a directory at 'path'."""
        if os.path.isdir(path):
            print('PATH: ', path, ' already exists')
        else:
            os.mkdir(path)
            print('PATH: ', path, ' created')


# =============================================================================
# 2dfdr aaorun wrapper
# =============================================================================
def aaorun_cleanup(log=False):
    """blah..."""
    process = subprocess.run('cleanup', shell=True, timeout=60, stdout=subprocess.PIPE, text=True)
    if log:
        logging.info('· [aaorun] Cleaning')


def aaorun_command(command, file, options=None, output=None,
                   idx_file='koala.idx',
                   aaorun='aaorun', wdir=None, timeout=900, log=False):
    """Create a tramline file from fibre flat.

    Link to 2dfdr repo:
    https://dev.aao.org.au/rds/2dfdr

    params
    ----
    file: (str) abs_path to file(s).
    output: (str) file name where the output will be stored.
    idx_file: (str) path to the idx file used for reducing data. Default is koala.idx
    aaorun: (str) path to the binary executable file
    """
    if wdir is None:
        wdir = os.path.dirname(file)

    if options is None:
        cmd_options = ' '.join(['-idxfile %s' % idx_file, '-wdir %s' % wdir])
    else:
        extra_options = [opt for opt in options]
        extra_options.append('-idxfile %s' % idx_file)
        extra_options.append('-wdir %s' % wdir)
        cmd_options = ' '.join(extra_options)
    # Combine all command arguments
    aaorun_cmd = ' '.join([aaorun, command, file, cmd_options])
    # Run the command
    if output is None:
        try:
            process = subprocess.run(aaorun_cmd, shell=True,
                                     timeout=timeout, stdout=subprocess.PIPE,
                                     text=True)
            aaorun_cleanup(log=log)
            if process.returncode != 0:
                if log:
                    logging.warning('[aaorun] · WARNING: Command \n {} \n FAILED! \n {}'.
                          format(aaorun_cmd, process.stderr))
                else:
                    print('[aaorun] · WARNING: Command \n {} \n FAILED! \n {}'.
                          format(aaorun_cmd, process.stderr))
                return 1
            else:
                return 0
        except subprocess.TimeoutExpired:
            if log:
                logging.warning('[aaorun] · WARNING: Command \n {} \n  Ran too long'.format(aaorun_cmd))
            else:
                print('[aaorun] · WARNING: Command \n {} \n  Ran too long'.format(aaorun_cmd))
            aaorun_cleanup(log=log)
            return 1
    else:
        with open(output, 'w') as outfile:
            try:
                process = subprocess.run(aaorun_cmd, shell=True,
                                         timeout=timeout, stdout=outfile,
                                         text=True)
                aaorun_cleanup(log=log)
                if process.returncode != 0:
                    if log:
                        logging.warning('[aaorun] · WARNING: Command \n {} \n FAILED! \n {}'.
                                        format(aaorun_cmd, process.stderr))
                    else:
                        print('[aaorun] · WARNING: Command \n {} \n FAILED! \n {}'.
                              format(aaorun_cmd, process.stderr))
                    return 1
                else:
                    return 0
            except subprocess.TimeoutExpired:
                if log:
                    logging.warning('[aaorun] · WARNING: Command \n {} \n  Ran too long'.format(aaorun_cmd))
                else:
                    print('[aaorun] · WARNING: Command \n {} \n  Ran too long'.format(aaorun_cmd))
                aaorun_cleanup(log=log)
                return 1


if __name__ == '__main__':
    aaorun_command('reduce_bias',
                   file='/media/pablo/toshiba-pab/test/20201021/ccd_1/21oct10001.fits',
                   aaorun='aaorun')

# Mr Krtxo
