import os
import yaml
import logging
import datetime
from typing import Dict, List, Optional
from pathlib import Path

import numpy as np

from py2dfdr.py2dfdr import aaorun_command
from py2dfdr import reduction_QC as QC
from py2dfdr import koala_cheatsheet as kcs
from py2dfdr import verbose


class ReduceObsRun:
    """Reduce KOALA observing runs.

    This class performs the main reduction steps and QC checks to reduce
    a KOALA observing run.

    Attributes
    ----------
    obs_run_info : dict
        Parsed contents of the OR 'obs_run_info.yml'.
    obs_run_path : str
        Root path to the observing run.
    ccds : list[str]
        CCDs to reduce (e.g., ['ccd_1', 'ccd_2']).
    dark_idx_file : Optional[str]
        .idx file for dark reductions.
    lflat_idx_file : Optional[str]
        .idx file for long-slit flat reductions.
    fibreflat_idx_file : Optional[str]
        .idx file for fibre-flat (incl. TLM) reductions.
    arc_idx_file : Optional[str]
        .idx file for arc reductions.
    object_idx_file : Optional[str]
        .idx file for science reductions.
    nights : list[str]
        List of observing nights in the OR (excludes calibration blocks).
    master_bias, master_darks, master_lflats, master_tlm, master_arcs,
    master_fibreflats : dict
        Paths to master calibration products, structured by CCD/night/grating.
    sat_fraction : float
        Max fraction of saturated pixels allowed (>= triggers rejection).
    sat_level : float
        Saturation level in counts.
    reject_names : list[str]
        Name keywords that trigger rejection (case-insensitive).
    """

    def __init__(self, obs_run_path: str, verb: bool = True, **kwargs):
        """
        Parameters
        ----------
        obs_run_path : str
            Path to the OR directory containing 'obs_run_info.yml'.
        verb : bool, optional
            If True, echo logs to stdout in addition to file logging.
        kwargs :
            sat_frac : float
            sat_level : float
            reject_names : List[str]
            ccds : List[str]              (required)
            dark_idx : Optional[str]
            lflat_idx : Optional[str]
            fibreflat_idx : Optional[str]
            arcs_idx : Optional[str]
            object_idx : Optional[str]
            master_dark_exptime : float   (target exptime when choosing a master)
            skip_nights : List[str]
            night_to_remove : List[str]   (legacy alias)
        """
        # Root path
        self.obs_run_path = os.fspath(obs_run_path)

        # Initialise logging (needs obs_run_path to exist)
        self.initialise_logger(verb=verb)

        # Data rejection thresholds & rules
        self.sat_level: float = float(kwargs.get('sat_level', 65500.0))
        self.sat_fraction: float = float(kwargs.get('sat_frac', 0.2))
        self.reject_names: List[str] = list(kwargs.get('reject_names', ['FOCUS']))

        # CCDs to process
        self.ccds: Optional[List[str]] = kwargs.get('ccds')
        if not self.ccds:
            logging.error('[OBSRUN] ERROR: CCDs not provided. Expected like ["ccd_1", "ccd_2"].')
            raise NameError('CCDs not provided (e.g. ["ccd_1", "ccd_2"]).')
        logging.info('CCDs to reduce: %s', ', '.join(self.ccds))

        # 2dfdr parameter files
        logging.info('[OBSRUN] Setting configuration files for 2dfdr')
        self.dark_idx_file: Optional[str] = kwargs.get('dark_idx')
        logging.info('[OBSRUN] DARK idx: %s', self.dark_idx_file)
        self.lflat_idx_file: Optional[str] = kwargs.get('lflat_idx')
        logging.info('[OBSRUN] LFLAT idx: %s', self.lflat_idx_file)
        self.fibreflat_idx_file: Optional[str] = kwargs.get('fibreflat_idx')
        logging.info('[OBSRUN] FFLAT idx: %s', self.fibreflat_idx_file)
        self.arc_idx_file: Optional[str] = kwargs.get('arcs_idx')
        logging.info('[OBSRUN] ARCS idx: %s', self.arc_idx_file)
        self.object_idx_file: Optional[str] = kwargs.get('object_idx')
        logging.info('[OBSRUN] OBJECT idx: %s', self.object_idx_file)

        # Targets for master selection
        self.master_dark_exptime: float = float(kwargs.get('master_dark_exptime', 1800))

        # Master containers
        self.obs_run_info: Optional[dict] = {}
        self.nights: Optional[List[str]] = None
        self.master_bias: Optional[dict] = None
        self.master_darks: Optional[Dict[str, str]] = None
        self.master_lflats: Optional[dict] = None
        self.master_tlm: Optional[dict] = None
        self.master_arcs: Optional[dict] = None
        self.master_fibreflats: Optional[dict] = None

        # Load OR descriptor
        self.load_obs_run_yml(**kwargs)
        # Grating to process
        self.gratings: Optional[List[str]] = kwargs.get("gratings",
                                                        self.get_gratings())
        # Optionally skip nights
        skip = kwargs.get("skip_nights", kwargs.get("night_to_remove", [])) or []
        if skip:
            logging.info("Skipping the following nights (if present): %s", ', '.join(skip))
            self.remove_nights(skip)
        # Optionally skip gratings
        skip = kwargs.get("skip_gratings", kwargs.get("gratings_to_remove", [])) or []
        if skip:
            logging.info("Skipping the following gratings (if present): %s", ', '.join(skip))
            self.remove_gratings(skip)
        
        self.skip_exist = kwargs.get("skip_exist", False)

    # ----------------------------------------------------------------------------------
    # Logging / IO helpers
    # ----------------------------------------------------------------------------------
    def initialise_logger(self, verb: bool) -> None:
        """Initialise file logger (and optional console echo)."""
        log_path = Path(self.obs_run_path) / 'OR_reduction.log'
        log_path.parent.mkdir(parents=True, exist_ok=True)

        logging.basicConfig(
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            filename=os.fspath(log_path),
            level=logging.INFO
        )
        if verb:
            logging.getLogger().addHandler(logging.StreamHandler())

        verbose.log_header(f'[OBSRUN] Initialising OR reduction at:\n   {self.obs_run_path}\n')
        logging.info('Session start: %s', datetime.datetime.now().strftime("%c"))

    def _read_recommended_path(self, file_path: str) -> Optional[str]:
        """Read first non-comment, non-empty line from a 'RECOMMENDED_*' file."""
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    s = line.strip()
                    if not s or s.startswith('#'):
                        continue
                    return s
        except Exception as e:
            logging.warning('Failed reading recommended file %s: %s', file_path, e)
        return None

    # ----------------------------------------------------------------------------------
    # OR descriptor
    # ----------------------------------------------------------------------------------
    def load_obs_run_yml(self, **kwargs) -> None:
        """Load the Observing Run YAML file and populate night list."""
        yml_path = Path(self.obs_run_path) / "obs_run_info.yml"
        logging.info('[OBSRUN] Loading OR descriptor: %s', yml_path)
        try:
            with open(yml_path, "r") as stream:
                self.obs_run_info = yaml.safe_load(stream) or {}
        except Exception as exc:
            logging.error('[OBSRUN] ERROR: Unable to load %s', yml_path)
            logging.error('%s', exc)
            raise

        # Nights present (exclude calibration blocks)
        raw_nights = list(self.obs_run_info.keys())
        calib_blocks = {'bias', 'darks', 'fflats'}
        self.nights = [n for n in raw_nights if n not in calib_blocks]

    # ----------------------------------------------------------------------------------
    # Master presence checks
    # ----------------------------------------------------------------------------------
    def check_masters(self, names: Optional[List[str]] = None) -> None:
        """If a master dict is missing, call its getter."""
        masters = {
            'bias': self.master_bias,
            'darks': self.master_darks,
            'lflats': self.master_lflats,
            'tlm': self.master_tlm,
            'arcs': self.master_arcs,
            'fflats': self.master_fibreflats
        }
        getters = {
            'bias': self.get_master_bias,
            'darks': self.get_master_darks,
            'lflats': self.get_master_lflats,
            'tlm': self.get_master_tlm,
            'arcs': self.get_master_arcs,
            'fflats': self.get_master_fibreflats
        }
        if not names:
            return
        for name in names:
            if masters.get(name) is None:
                getters[name]()

    # ----------------------------------------------------------------------------------
    # Rejection rules
    # ----------------------------------------------------------------------------------
    def reject_saturated_image(self, path: str) -> bool:
        """Reject object based on fraction of saturated pixels."""
        sat_frac = QC.check_saturated(path, sat_level=self.sat_level, log=True)
        return bool(sat_frac >= self.sat_fraction)

    def reject_from_name(self, name: str) -> bool:
        """Reject object based on header name."""
        lname = name.lower()
        return any(r.lower() in lname for r in self.reject_names)

    @staticmethod
    def reject_binned(path: str) -> bool:
        """Reject data taken in a binned/windowed mode."""
        val = QC.get_keyword(path, keyword='WINDOW', hdu_index=0)
        if val is None:
            return False
        return 'BIN' in str(val)

    def select_nights(self, nights: List[str]) -> None:
        """Manually set nights to be reduced (overrides autodetected)."""
        logging.info("Manually setting nights to reduce: %s", ', '.join(nights))
        self.nights = list(nights)

    def remove_nights(self, nights: List[str]) -> None:
        """Remove given nights from the reduction list (if present)."""
        if not self.nights:
            return
        keep = set(self.nights) - set(nights)
        removed = set(self.nights) & set(nights)
        self.nights = [n for n in self.nights if n in keep]
        if removed:
            for n in sorted(removed):
                logging.info("Skipping night: %s", n)

    def get_gratings(self):
        """Retrieve all gratings included in the reduction list."""
        logging.info("Retrieving all gratings")
        if not self.nights:
            return
        gratings = []
        for night in self.nights or []:
            for ccd in self.ccds:
                if ccd not in self.obs_run_info.get(night, {}):
                    continue
                for grating in self.obs_run_info[night][ccd].keys():
                    if grating not in gratings:
                        gratings.append(grating)
        logging.info("List of selected gratings: " + ", ".join(gratings))
        return gratings

    def remove_gratings(self, gratings: List[str]) -> None:
        if not self.gratings:
            return
        keep = set(self.gratings) - set(gratings)
        removed = set(self.gratings) & set(gratings)
        self.gratings = [g for g in self.gratings if g in keep]
        if removed:
            for g in sorted(removed):
                logging.info("Skipping grating: %s", g)

    # ----------------------------------------------------------------------------------
    # Reduction steps
    # ----------------------------------------------------------------------------------
    def reduce_bias(self) -> None:
        """(Placeholder) Bias handling if needed in the future."""
        verbose.log_header('bias')
        if not self.obs_run_info or 'bias' not in self.obs_run_info:
            logging.info('No bias block found.')
            return
        for ccd in self.ccds:
            if ccd not in self.obs_run_info['bias']:
                continue
            for _, bias_file in self.obs_run_info['bias'][ccd].items():
                logging.info('BIAS: %s', bias_file.get('PATH'))

    def get_master_bias(self) -> None:
        """(Placeholder) Locate master bias if used by your workflow."""
        self.master_bias = {}

    def reduce_darks(self, timeout: int = 900) -> None:
        """Reduce dark frames and build master darks."""
        verbose.log_header('Reducing darks')
        if not self.dark_idx_file:
            verbose.missing_idx('dark_idx')

        for ccd in self.ccds:
            if 'darks' not in self.obs_run_info or ccd not in self.obs_run_info['darks']:
                continue
            for exptime in self.obs_run_info['darks'][ccd].keys():
                all_darks = []
                for _, dark_file in self.obs_run_info['darks'][ccd][exptime].items():
                    path_to_dark = os.path.join(
                        self.obs_run_path, 'darks', ccd, exptime, dark_file['PATH'])
                    logging.info('[OBSRUN] [%s] [%s] DARK: %s', ccd, exptime, path_to_dark)
                    if not QC.check_exists(path_to_dark):
                        logging.warning('File does not exist: %s', path_to_dark)
                        continue
                    aaorun_command('reduce_dark', path_to_dark,
                                   idx_file=self.dark_idx_file,
                                   output=path_to_dark.replace('.fits', '_log.txt'),
                                   timeout=timeout)
                    all_darks.append(path_to_dark.replace('.fits', 'red.fits'))

                if not all_darks:
                    logging.warning('[OBSRUN] No darks to combine for [%s] [%s]', ccd, exptime)
                    continue

                logging.info('[OBSRUN] Computing MASTERDARK [%s] [%s]', ccd, exptime)
                darks_to_combine = ' '.join(all_darks)
                masterdark_name = os.path.join(
                    self.obs_run_path, 'darks', ccd, exptime, f'DARKcombined_{exptime}.fits')
                logmaster = masterdark_name.replace('.fits', '_log.txt')

                aaorun_command(command='combine_image',
                               file=f'"{darks_to_combine}"',
                               idx_file=self.dark_idx_file,
                               options=[f'-COMBINEDFILE {masterdark_name}'],
                               output=logmaster,
                               wdir=os.path.dirname(masterdark_name),
                               log=True)
                logging.info('[OBSRUN] MASTERDARK: %s', masterdark_name)
                QC.check_image(masterdark_name,
                               save_dir=os.path.join(self.obs_run_path, 'darks', ccd, exptime, 'masterdark.png'))

        # Select best master per CCD by nearest exptime (abs diff)
        self.master_darks = {}
        for ccd in self.ccds:
            if 'darks' not in self.obs_run_info or ccd not in self.obs_run_info['darks']:
                continue
            exptimes = list(self.obs_run_info['darks'][ccd].keys())
            if not exptimes:
                continue
            arr = np.array(exptimes, dtype=float)
            bestexp = int(np.argmin(np.abs(arr - self.master_dark_exptime)))
            optimal_exptime = exptimes[bestexp]
            self.master_darks[ccd] = os.path.join(
                self.obs_run_path, 'darks', ccd, optimal_exptime, f'DARKcombined_{optimal_exptime}.fits'
            )

    def get_master_darks(self, exptime: Optional[float] = None) -> None:
        """Locate existing master darks on disk (per CCD)."""
        self.master_darks = {}
        verbose.log_header('[OBSRUN] Searching master dark files')
        for ccd in self.ccds:
            if 'darks' not in self.obs_run_info or ccd not in self.obs_run_info['darks']:
                continue
            exptimes = list(self.obs_run_info['darks'][ccd].keys())
            if not exptimes:
                continue
            target = self.master_dark_exptime if exptime is None else float(exptime)
            arr = np.array(exptimes, dtype=float)
            bestexp = int(np.argmin(np.abs(arr - target)))
            optimal_exptime = exptimes[bestexp]
            path_to_master = os.path.join(
                self.obs_run_path, 'darks', ccd, optimal_exptime, f'DARKcombined_{optimal_exptime}.fits')
            if os.path.isfile(path_to_master):
                logging.info('[OBSRUN] [%s] MASTERDARK: %s', ccd, path_to_master)
                self.master_darks[ccd] = path_to_master
            else:
                verbose.missing_master(path_to_master)

    def reduce_lflats(self, timeout: int = 900) -> None:
        """Reduce long-slit flats and select a recommended master per CCD+grating."""
        verbose.log_header('Reducing long-slit (detector) flats')
        if not self.lflat_idx_file:
            verbose.missing_idx('lflat_idx')

        self.master_lflats = {}
        for ccd in self.ccds:
            if 'fflats' not in self.obs_run_info or ccd not in self.obs_run_info['fflats']:
                continue
            gratings = set(self.obs_run_info['fflats'][ccd]) & set(self.gratings)
            self.master_lflats[ccd] = {}
            for grating in gratings:
                self.master_lflats[ccd][grating] = {}
                for name, fflat_file in self.obs_run_info['fflats'][ccd][grating].items():
                    fflat_exptime = fflat_file['EXPTIME']
                    self.master_lflats[ccd][grating].setdefault(fflat_exptime, [])

                    path_to_fflat = os.path.join(self.obs_run_path, 'fflats', ccd, grating, fflat_file['PATH'])
                    logging.info('[OBSRUN] [%s] [%s] LFLAT: %s', ccd, grating, path_to_fflat)

                    if not QC.check_exists(path_to_fflat):
                        logging.warning('File does not exist: %s', path_to_fflat)
                        continue

                    aaorun_command('reduce_lflat', path_to_fflat,
                                   idx_file=self.lflat_idx_file,
                                   output=path_to_fflat.replace('.fits', '_log.txt'),
                                   log=True, timeout=timeout)

                    self.master_lflats[ccd][grating][fflat_exptime].append(
                        path_to_fflat.replace('.fits', 'red.fits')
                    )

                # Combine per exposure time
                for exptime, files in list(self.master_lflats[ccd][grating].items()):
                    if len(files) <= 1:
                        # Keep the single file as the "master" for that exptime
                        if files:
                            self.master_lflats[ccd][grating][exptime] = files[0]
                            logging.warning('[OBSRUN] Only one LFLAT for [%s] [%s] [%s]; skipping combine.',
                                            ccd, grating, exptime)
                        continue

                    verbose.log_header(f'[OBSRUN] Computing MASTERLFLAT [{ccd}] [{grating}] [{exptime}]')
                    fflats_to_combine = ' '.join(files)
                    masterfflat_name = os.path.join(
                        self.obs_run_path, 'fflats', ccd, grating, f'LFLATcombined_{grating}_{exptime}.fits')
                    logmaster = masterfflat_name.replace('.fits', '_log.txt')

                    aaorun_command(command='combine_image',
                                   file=f'"{fflats_to_combine}"',
                                   options=[f'-COMBINEDFILE {masterfflat_name}'],
                                   output=logmaster,
                                   wdir=os.path.dirname(masterfflat_name),
                                   log=True)

                    self.master_lflats[ccd][grating][exptime] = masterfflat_name
                    logging.info('[OBSRUN] MASTERLFLAT: %s', masterfflat_name)
                    QC.check_image(masterfflat_name, save_dir=os.path.join(
                        self.obs_run_path, 'fflats', ccd, grating, f'masterlflat_{exptime}.png'))

                # Pick recommended by closest to KOALA cheatsheet time
                try:
                    rec_time = kcs.lflat_time[ccd][grating.split('_')[0]]
                    time_keys = list(self.master_lflats[ccd][grating].keys())
                    times = np.array(time_keys, dtype=float)
                    best = int(np.argmin(np.abs(times - rec_time)))
                    selected = self.master_lflats[ccd][grating][time_keys[best]]
                    self.master_lflats[ccd][grating] = selected
                    rec_file = os.path.join(os.path.dirname(selected), 'RECOMMENDED_FLAT')
                    with open(rec_file, 'w') as f:
                        f.write(str(selected))
                    logging.info('[OBSRUN] [%s] [%s] LFLAT recommended: %s (target %.1f s)',
                                 ccd, grating, selected, rec_time)
                except Exception as e:
                    logging.warning('[OBSRUN] Failed to select recommended LFLAT for [%s][%s]: %s',
                                    ccd, grating, e)

    def get_master_lflats(self) -> None:
        """Load selected (recommended) master LFLAT per CCD+grating."""
        self.master_lflats = {}
        logging.info('[OBSRUN] Searching master LFLAT files')
        for ccd in self.ccds:
            if 'fflats' not in self.obs_run_info or ccd not in self.obs_run_info['fflats']:
                continue
            gratings = set(self.obs_run_info['fflats'][ccd]) & set(self.gratings)
            self.master_lflats[ccd] = {}
            for grating in gratings:
                rec_path = os.path.join(self.obs_run_path, 'fflats', ccd, grating, 'RECOMMENDED_FLAT')
                if not os.path.isfile(rec_path):
                    verbose.NoFileError(rec_path)
                    continue
                path_to_master = self._read_recommended_path(rec_path)
                if path_to_master and os.path.isfile(path_to_master):
                    logging.info('[OBSRUN] [%s] [%s] MASTERLFLAT: %s', ccd, grating, path_to_master)
                    self.master_lflats[ccd][grating] = path_to_master
                else:
                    verbose.missing_master(path_to_master or rec_path)

    def extract_tramlines(self, timeout: int = 900) -> None:
        """Extract tramline maps from fibre flats, select recommended per night/CCD/grating."""
        verbose.log_header('Starting tramline extraction from fibre flats')
        if not self.fibreflat_idx_file:
            verbose.missing_idx('fflat_idx')

        self.master_tlm = {}
        for night in self.nights or []:
            self.master_tlm[night] = {}
            for ccd in self.ccds:
                self.master_tlm[night][ccd] = {}
                if ccd not in self.obs_run_info.get(night, {}):
                    continue
                for grating in set(self.obs_run_info[night][ccd]) & set(self.gratings):
                    exptimes: List[float] = []
                    names: List[str] = []
                    for name, tram_file in self.obs_run_info[night][ccd][grating]['fibreflat'].items():
                        path_to_fibreflat = os.path.join(
                            self.obs_run_path, night, ccd, grating, 'fibreflat', tram_file['PATH'])
                        logging.info('[OBSRUN] [%s] [%s] [%s] TRAM: %s', night, ccd, grating, path_to_fibreflat)

                        if not QC.check_exists(path_to_fibreflat):
                            logging.warning('File does not exist: %s', path_to_fibreflat)
                            continue
                        if self.reject_saturated_image(path_to_fibreflat):
                            logging.info('[OBSRUN] [%s] [%s] [%s] TRAM: SATURATED -> skip.', night, ccd, grating)
                            continue

                        ok = aaorun_command('make_tlm', path_to_fibreflat,
                                            idx_file=self.fibreflat_idx_file,
                                            output=path_to_fibreflat.replace('.fits', '_log.txt'),
                                            log=True, timeout=timeout)
                        if ok == 0:
                            bad_tlm = QC.check_tramline(path_to_fibreflat.replace('.fits', 'tlm.fits'))
                        else:
                            bad_tlm = True

                        if not bad_tlm:
                            exptimes.append(float(tram_file['EXPTIME']))
                            names.append(name)

                    try:
                        rec_time = kcs.lflat_time[ccd][grating.split('_')[0]]
                    except Exception:
                        rec_time = None

                    if exptimes:
                        idx = int(np.argmin(np.abs(np.array(exptimes) - (rec_time or np.median(exptimes)))))
                        best_fibreflat = os.path.join(
                            self.obs_run_path, night, ccd, grating, 'fibreflat',
                            self.obs_run_info[night][ccd][grating]['fibreflat'][names[idx]]['PATH']
                        ).replace('.fits', 'tlm.fits')

                        logging.info('[OBSRUN] [%s] [%s] [%s] BEST TRAM: %s', night, ccd, grating, best_fibreflat)
                        self.master_tlm[night][ccd][grating] = best_fibreflat
                        rec_file = os.path.join(os.path.dirname(best_fibreflat), 'RECOMMENDED_TRAM')
                        with open(rec_file, 'w') as f:
                            f.write(best_fibreflat)
                    else:
                        logging.warning('[OBSRUN] [%s] [%s] [%s] WARNING: NO TRAMLINE AVAILABLE.', night, ccd, grating)

    def get_master_tlm(self) -> None:
        """Locate tramline maps from RECOMMENDED_TRAM files."""
        self.master_tlm = {}
        verbose.log_header('[OBSRUN] Searching TRAMLINE MAPS')
        for night in self.nights or []:
            self.master_tlm[night] = {}
            for ccd in self.ccds:
                self.master_tlm[night][ccd] = {}
                for grating in set(self.obs_run_info[night][ccd]) & set(self.gratings):
                    rec_path = os.path.join(self.obs_run_path, night, ccd, grating, 'fibreflat', 'RECOMMENDED_TRAM')
                    if not os.path.isfile(rec_path):
                        verbose.NoFileError(rec_path)
                        continue
                    path_to_tram = self._read_recommended_path(rec_path)
                    if path_to_tram and os.path.isfile(path_to_tram):
                        logging.info('[OBSRUN] [%s] [%s] [%s] TRAM: %s', night, ccd, grating, path_to_tram)
                        self.master_tlm[night][ccd][grating] = path_to_tram
                    else:
                        verbose.missing_master(path_to_tram or rec_path)

    def reduce_arcs(self, timeout: int = 900) -> None:
        """Reduce calibration arcs and pick recommended per night/CCD/grating."""
        verbose.log_header('[OBSRUN] Reducing calibration arcs')
        if not self.arc_idx_file:
            verbose.missing_idx('arcs_idx')

        self.master_arcs = {}
        self.check_masters(names=['darks', 'lflats', 'tlm'])
        for night in self.nights or []:
            self.master_arcs[night] = {}
            for ccd in self.ccds:
                self.master_arcs[night][ccd] = {}
                if ccd not in self.obs_run_info.get(night, {}):
                    continue
                for grating in set(self.obs_run_info[night][ccd]) & set(self.gratings):
                    exptimes: List[float] = []
                    names: List[str] = []
                    arcnames: List[str] = []

                    for name, arc_file in self.obs_run_info[night][ccd][grating]['arcs'].items():
                        arc_name = arc_file['ARCNAME']
                        arc_exptime = float(arc_file['EXPTIME'])
                        path_to_arc = os.path.join(self.obs_run_path, night, ccd, grating, 'arcs', arc_file['PATH'])
                        logging.info('[OBSRUN] [%s] [%s] [%s] ARC %s (%.1f s): %s',
                                     night, ccd, grating, arc_name, arc_exptime, path_to_arc)

                        if not QC.check_exists(path_to_arc):
                            logging.warning('File does not exist: %s', path_to_arc)
                            continue
                        if self.reject_saturated_image(path_to_arc):
                            logging.info('[OBSRUN] [%s] [%s] [%s] ARC SATURATED -> skip.', night, ccd, grating)
                            continue

                        ok = aaorun_command(
                            'reduce_arc', path_to_arc,
                            idx_file=self.arc_idx_file,
                            options=[
                                f'-DARK_FILENAME {self.master_darks.get(ccd, "")}',
                                f'-TLMAP_FILENAME {self.master_tlm.get(night, {}).get(ccd, {}).get(grating, "")}'
                            ],
                            output=path_to_arc.replace('.fits', '_log.txt'),
                            log=True, timeout=timeout
                        )
                        if ok == 0:
                            QC.check_image(path_to_arc.replace('.fits', 'red.fits'),
                                           title=f"{Path(path_to_arc).name}\narc: {arc_name}-{arc_exptime:.0f} sec",
                                           save_dir=path_to_arc.replace('.fits', '_qc.png'))
                            names.append(name)
                            exptimes.append(arc_exptime)
                            arcnames.append(arc_name)

                    # Select recommended by name+closest exptime
                    try:
                        rec_arc_name, rec_exptime = kcs.arc_time[ccd][grating.split('_')[0]]
                    except Exception:
                        rec_arc_name, rec_exptime = None, None

                    if names:
                        selected_mask = np.array([rec_arc_name in a if rec_arc_name else True for a in arcnames])
                        expt_arr = np.array(exptimes, dtype=float)
                        expt_arr[~selected_mask] = np.nan
                        if np.isfinite(expt_arr).any():
                            best_local = int(np.nanargmin(np.abs(expt_arr[selected_mask] - (rec_exptime or np.nanmedian(expt_arr)))))
                            selected_indices = np.where(selected_mask)[0]
                            best_idx = int(selected_indices[best_local])
                        else:
                            # Fallback: choose global closest to median exptime
                            best_idx = int(np.argmin(np.abs(expt_arr - np.nanmedian(expt_arr))))

                        best_arc = os.path.join(
                            self.obs_run_path, night, ccd, grating, 'arcs',
                            self.obs_run_info[night][ccd][grating]['arcs'][names[best_idx]]['PATH']
                        ).replace('.fits', 'red.fits')

                        logging.info('[OBSRUN] [%s] [%s] [%s] BEST ARC: %s', night, ccd, grating, best_arc)
                        self.master_arcs[night][ccd][grating] = best_arc
                        rec_file = os.path.join(os.path.dirname(best_arc), 'RECOMMENDED_ARC')
                        with open(rec_file, 'w') as f:
                            f.write(best_arc)
                    else:
                        logging.warning('[OBSRUN] [%s] [%s] [%s] WARNING: NO ARC SELECTED!', night, ccd, grating)

    def get_master_arcs(self) -> None:
        """Locate reduced arcs from RECOMMENDED_ARC files."""
        self.master_arcs = {}
        verbose.log_header('[OBSRUN] Searching ARCS')
        for night in self.nights or []:
            self.master_arcs[night] = {}
            for ccd in self.ccds:
                self.master_arcs[night][ccd] = {}
                for grating in set(self.obs_run_info[night][ccd]) & set(self.gratings):
                    rec_path = os.path.join(self.obs_run_path, night, ccd, grating, 'arcs', 'RECOMMENDED_ARC')
                    if not os.path.isfile(rec_path):
                        verbose.NoFileError(rec_path)
                        continue
                    path_to_arc = self._read_recommended_path(rec_path)
                    if path_to_arc and os.path.isfile(path_to_arc):
                        logging.info('[OBSRUN] [%s] [%s] [%s] ARC: %s', night, ccd, grating, path_to_arc)
                        self.master_arcs[night][ccd][grating] = path_to_arc
                    else:
                        verbose.missing_master(path_to_arc or rec_path)

    def reduce_fflats(self, timeout: int = 900) -> None:
        """Reduce fibre flats (using selected TLM, DARK, ARC) and store masters."""
        verbose.log_header('[OBSRUN] Reducing FIBRE FLATS')
        if not self.fibreflat_idx_file:
            verbose.missing_idx('fflat_idx')

        self.check_masters(names=['darks', 'lflats', 'tlm', 'arcs'])

        # Seed master_fibreflats structure based on master_tlm keys, replacing tlm->red
        self.master_fibreflats = {}
        for night, d_ccd in (self.master_tlm or {}).items():
            self.master_fibreflats[night] = {}
            for ccd, d_gr in d_ccd.items():
                self.master_fibreflats[night][ccd] = {}
                for grating, tlm_path in d_gr.items():
                    self.master_fibreflats[night][ccd][grating] = tlm_path.replace('tlm.fits', 'red.fits')

        for night in self.nights or []:
            for ccd in self.ccds:
                if ccd not in self.obs_run_info.get(night, {}):
                    continue
                for grating in set(self.obs_run_info[night][ccd]) & set(self.gratings):
                    for name, fflat in self.obs_run_info[night][ccd][grating]['fibreflat'].items():
                        path_to_fflat = os.path.join(self.obs_run_path, night, ccd, grating, 'fibreflat', fflat['PATH'])
                        exptime = fflat['EXPTIME']
                        logging.info('[OBSRUN] [%s] [%s] [%s] FFLAT %s s: %s', night, ccd, grating, exptime, path_to_fflat)

                        if not QC.check_exists(path_to_fflat):
                            logging.warning('File does not exist: %s', path_to_fflat)
                            continue
                        if self.reject_saturated_image(path_to_fflat):
                            logging.info('[OBSRUN] [%s] [%s] [%s] FFLAT SATURATED -> skip.', night, ccd, grating)
                            continue

                        ok = aaorun_command(
                            'reduce_fflat', path_to_fflat,
                            idx_file=self.fibreflat_idx_file,
                            options=[
                                f'-DARK_FILENAME {self.master_darks.get(ccd, "")}',
                                f'-TLMAP_FILENAME {self.master_tlm.get(night, {}).get(ccd, {}).get(grating, "")}',
                                f'-WAVEL_FILENAME {self.master_arcs.get(night, {}).get(ccd, {}).get(grating, "")}'
                            ],
                            output=path_to_fflat.replace('.fits', '_log.txt'),
                            log=True, timeout=timeout
                        )
                        if ok == 0:
                            out = path_to_fflat.replace('.fits', 'red.fits')
                            QC.check_image(out, save_dir=os.path.join(os.path.dirname(path_to_fflat), f'{name}.png'))
                            # Keep the one that corresponds to the chosen TLM (narrow heuristic)
                            self.master_fibreflats.setdefault(night, {}).setdefault(ccd, {})[grating] = out
                        else:
                            logging.warning('[OBSRUN] [%s] [%s] [%s] FFLAT reduction failed.', night, ccd, grating)

    def get_master_fibreflats(self) -> None:
        """Locate reduced fibre flats based on TRAM recommendations."""
        self.master_fibreflats = {}
        verbose.log_header('[OBSRUN] Searching fibre flats')
        for night in self.nights or []:
            self.master_fibreflats[night] = {}
            for ccd in self.ccds:
                self.master_fibreflats[night][ccd] = {}
                for grating in set(self.obs_run_info[night][ccd]) & set(self.gratings):
                    rec_tram = os.path.join(self.obs_run_path, night, ccd, grating, 'fibreflat', 'RECOMMENDED_TRAM')
                    if not os.path.isfile(rec_tram):
                        verbose.NoFileError(rec_tram)
                        continue
                    path_to_tram = self._read_recommended_path(rec_tram)
                    if not path_to_tram:
                        verbose.missing_master(rec_tram)
                        continue
                    path_to_fflat = path_to_tram.replace('tlm.fits', 'red.fits')
                    if os.path.isfile(path_to_fflat):
                        logging.info('[OBSRUN] [%s] [%s] [%s] FFLAT: %s', night, ccd, grating, path_to_fflat)
                        self.master_fibreflats[night][ccd][grating] = path_to_fflat
                    else:
                        verbose.missing_master(path_to_fflat)

    def reduce_object(self, timeout: int = 900) -> None:
        """Reduce science exposures using selected masters and write REDUCTION_FLAGS.yml."""
        verbose.log_header('[OBSRUN] Reducing science objects')
        if not self.object_idx_file:
            verbose.missing_idx('object_idx')

        self.check_masters(names=['darks', 'lflats', 'tlm', 'arcs', 'fflats'])

        for night in self.nights or []:
            for ccd in self.ccds:
                if ccd not in self.obs_run_info.get(night, {}):
                    continue
                for grating in set(self.obs_run_info[night][ccd]) & set(self.gratings):
                    flags: Dict[str, dict] = {}
                    for name, obj in self.obs_run_info[night][ccd][grating]['sci'].items():
                        obj_name = obj['NAME']
                        exptime = obj['EXPTIME']
                        flags[name] = {'NAME': obj_name, 'PATH': obj['PATH']}

                        path_to_obj = os.path.join(self.obs_run_path, night, ccd, grating, 'sci', obj['PATH'])
                        logging.info('[OBSRUN] [%s] [%s] [%s] OBJ %s (%s s): %s',
                                     night, ccd, grating, obj_name, exptime, path_to_obj)
                        if not QC.check_exists(path_to_obj):
                            logging.warning('File does not exist: %s', path_to_obj)
                            continue
                        elif self.skip_exist:
                            if self.is_reduced(path_to_obj):
                                logging.warning('File is already reduced: SKIPPING')
                                continue
                        if self.reject_from_name(obj_name):
                            flags[name]['FLAG'] = 'NAMEREJ'
                            logging.info('[OBSRUN] Name-rejected -> skip.')
                            continue
                        if self.reject_binned(path_to_obj):
                            flags[name]['FLAG'] = 'BINNING_CONFIG'
                            logging.info('[OBSRUN] Binning config -> skip.')
                            continue
                        if self.reject_saturated_image(path_to_obj):
                            flags[name]['FLAG'] = 'SATURATED'
                            logging.info('[OBSRUN] Saturated -> skip.')
                            continue

                        ok = aaorun_command(
                            'reduce_object', path_to_obj,
                            idx_file=self.object_idx_file,
                            options=[
                                f'-DARK_FILENAME {self.master_darks.get(ccd, "")}',
                                f'-FFLAT_FILENAME {self.master_fibreflats.get(night, {}).get(ccd, {}).get(grating, "")}',
                                f'-TLMAP_FILENAME {self.master_tlm.get(night, {}).get(ccd, {}).get(grating, "")}',
                                f'-WAVEL_FILENAME {self.master_arcs.get(night, {}).get(ccd, {}).get(grating, "")}'
                            ],
                            output=path_to_obj.replace('.fits', '_log.txt'),
                            log=True, timeout=timeout
                        )
                        if ok == 0:
                            flags[name]['FLAG'] = 'OK'
                            QC.check_image(path_to_obj.replace('.fits', 'red.fits'),
                                           save_dir=path_to_obj.replace('.fits', '_qc.png'),
                                           title=obj_name)
                        else:
                            flags[name]['FLAG'] = 'AAORUNFAIL'
                            logging.warning('[OBSRUN] Reduction failed for: %s', path_to_obj)

                    # Persist flags summary
                    out_flags = os.path.join(self.obs_run_path, night, ccd, grating, 'sci', 'REDUCTION_FLAGS.yml')
                    try:
                        with open(out_flags, 'w') as outfile:
                            yaml.dump(flags, outfile, default_flow_style=False)
                        logging.info('[OBSRUN] Wrote flags: %s', out_flags)
                    except Exception as e:
                        logging.warning('[OBSRUN] Failed writing flags file %s: %s', out_flags, e)

    def is_reduced(self, path_to_obj):
        if os.path.isfile(path_to_obj.replace(".fits", "red.fits")):
            return True
        else:
            return False

    # ----------------------------------------------------------------------------------
    # Extra
    # ----------------------------------------------------------------------------------
    def combine_science_data(self, keyword: Optional[str] = None) -> None:
        """Combine all files within each night sharing a common name (TBD)."""
        pass
