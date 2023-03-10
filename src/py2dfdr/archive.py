import os
import numpy as np
from glob import glob
from astropy.io import fits
from astropy.time import Time
import shutil
import yaml
import logging


class ArchiveObs(object):
    """This class allows the user to classify and extract all the information from a set of observing runs.

    Attributes
    ----------
    - nights: (dict) Collection of observing nights containing file information.
    - dates: (array) Observing nights julian dates.
    - observing_runs: (dict) Collection of observing runs with file information.
    - root_path: (str) Path to the directory containing all subdirectories, containing the individual observing nights.

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

    def __init__(self, root_path, verbose=True):
        self.verbose = verbose
        self.archive_print("Initialising archive process...")
        # Declare the main attributes
        self.nights = None
        self.dates = None
        self.observing_runs = None
        self.root_path = root_path

    def archive_print(self, msg):
        """Print a given message if verbose is set to True.
        Params
        ------
        - msg: (str)
        """
        if self.verbose:
            print("[Archiving] " + msg)

    def find_nights(self):
        """Find individual folders within the root directory.
        
        Each individual directory is expected to follow the same name convention
        yr~month~day (e.g. 20180310) and are then converted to Julian Dates.
        """
        # Find elements within the root directory
        self.archive_print("Looking for observing nights at $root: \n{}".format(self.root_path))
        nights_path = np.sort(glob(os.path.join(self.root_path, "*")))
        self.archive_print("{} elements found withing $root".format(len(nights_path)))
        # Remove any file that could be within the root folder
        self.archive_print("Filtering files")
        nights_path = list(filter(os.path.isdir, nights_path))
        # Check if there is any night
        if len(nights_path) == 0:
            raise NameError(
                "Root path:\n{}\n contains no nights (e.g. $root/20180310)".format(self.root_path))
        else:
            self.archive_print("{} nights found".format(len(nights_path)))
            nights = {}
            # Julian Dates
            self.archive_print("Converting to Julian Dates")
            self.dates = np.zeros(len(nights_path))
            for i, night_dir in enumerate(nights_path):
                folder = night_dir.split('/')[-1]
                julian_date = Time(
                    folder[:4] + '-' + folder[4:6] + '-'+folder[6:]).jd
                self.dates[i] = julian_date
                nights[folder] = {'JD': julian_date}
            self.nights = nights

    def get_night_info(self):
        """Obtain the content of each observing night."""
        self.archive_print("Obtaning the content information of each night")
        for night in self.nights.keys():
            night_path = os.path.join(self.root_path, night)
            ccd_folders = glob(os.path.join(night_path, '*'))
            ccd_folders = list(filter(os.path.isdir, ccd_folders))
            self.nights[night]['content'] = {}
            for ccd in ccd_folders:
                content = self.get_content(ccd)
                self.nights[night]['content'][ccd.split('/')[-1]] = content
                ## Fill all nights with the same ORI
                self.nights[night]['ori'] = 0

    def get_content(self, path):
        """
        Return a dictionary containing the description of each file within a given path folder.

        params:
        ------
        - path (str): path to the data

        returns:
        --------
         - files: (dict)
        """
        ####
        # TODO: This function is hardcoded and should be provided by the user
        # based on the instrument in use (e.g. KOALA, SAMI, etc...)
        ####
        # Change the working directory
        os.chdir(path)
        self.archive_print('Searching files at: {}'.format(os.getcwd()))
        # Call aaorun list for selecting the elements within the folder
        aaorun_list = os.system('aaorun list > content.txt')
        if aaorun_list != 0:
            raise NameError('ERROR: aaorun list command did not work!')
        # Read the ouput from aaorun and create the "files" dictionary
        with open('content.txt', 'r') as f:
            content = f.read()
            content = content.split('\n')
            self.archive_print('Header line --> ' +  content[0])
            files_prefix = content[0][
                content[0].find(': ')+2: content[0].find(' (')]
            self.archive_print(' > Files prefix *{}*'.format(files_prefix))
            content = content[1:-1]
            files = {}
            for file in range(len(content)):
                data = content[file].split(' ')
                file_number = data[0].split(files_prefix)[1].split('.fits')[0]
                files[file_number] = {'PATH': data[0], 'DESC': data[1], 'REDUCE': True}
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
        self.archive_print(' > Found {} files'.format(len(files)))
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        return files

    def filter_night_files(self):
        raise NotImplementedError("Not implemented :(")

    def create_ori(self):
        """Assignate an observing run ID (ORI) to every observing night.
        
        Based on the Julian Date of each night, consecutive nights
        will be grouped together forming an observing run (see doc).

        Example
        -------
        If the $root dir contains the following observing night of
        a given month:
        1, 2, 3, 4, 8, 10, 11

        Then, the observing run ID (ORI) of each night is
        1, 2, 3, 4 --> 1
        8 --> 2
        10, 11 --> 3
        """
        self.archive_print("Grouping nights into observing runs")
        if len(self.dates) > 1:
            diffdates = np.diff(self.dates)
            binnumber = np.arange(0, len(self.dates), 1)
            # Bin the observing nights
            for i in range(diffdates.size):
                if diffdates[i] == 1:
                    binnumber[i+1] = binnumber[i]
            # Rename the bin numbers
            bins = np.unique(binnumber)
            self.ori = np.arange(0, bins.size)
            mapping = dict(zip(bins, self.ori))
            self.night_run_id = binnumber
            self.nights_keys = np.array(list(self.nights.keys()))
            for i, night in enumerate(self.nights.keys()):
                self.nights[night]['ori'] = mapping[binnumber[i]]
            self.archive_print("{} observing runs".format(bins.size))
        else:
            self.archive_print("All files correspond to a single observing run [ORI=0]")

    def create_reduction_folder(self, output_path=None, hard=True):
        """Create a reduction directory for each observing run, based on the available nights.
        
        Params
        ------
        output_path: (str) Reduction diractory root path
        hard: (bool) if True, hard copies of each file will be created. Otherwise
        symbolic links will be used to copy the files."""
        # TODO: The expected CCDs are hardcoded.
        self.archive_print("Creating reduction directory")
        # Create root dir for reduction
        if output_path is None:
            output_path = self.root_path
        self.archive_print('Reduction dir:\n '.format(output_path))
        if hard:
            copier = shutil.copy
            self.archive_print("Files will be hard-copied")
        else:
            copier = os.symlink
            self.archive_print("Files will be trated as symbolic links")

        reduce_path = os.path.join(output_path)
        self.make_dir(reduce_path)
        # Create observing runs folders
        self.observing_runs = {}
        for obs_run in self.ori:
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
        """Check whether a directory exists and create it otherwise.
        
        Params
        ------
        path: (str) path to the directory
        """
        if os.path.isdir(path):
            self.archive_print('PATH: ' +  path + ' already exists')
        else:
            os.mkdir(path)
            self.archive_print('PATH: ' + path + ' created')


