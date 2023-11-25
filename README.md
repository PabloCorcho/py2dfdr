# py2dfdr
2dfdr python wrapper for reducing AAT data.

[Link to AAT webpage](https://aat.anu.edu.au/)
[Link to 2dfdr webpage](https://aat.anu.edu.au/science/software/2dfdr)
[2dfdr repository](https://dev.aao.org.au/rds/2dfdr)

## Main scope of py2dfdr

This package aims to facilitate the data reduction process using the software 2dfdr following an object-oriented scheme that easily allows the data reduction of large amounts of data. The original purpose of this package is to reduce [KOALA](https://aat.anu.edu.au/science/instruments/current/koala/overview) data as part of the [HI-KIDS]() galaxy survey. Thus, it might not be suitable for reducing data from other AAT instruments (i.e., 2df, HERMES, SAMI, Hector, etc...)

Contact: (Pablo Corcho-Caballero) corchocaballeropablo@gmail.com / p.corcho.caballero@rug.nl

## Create observing runs data directory
Before starting the data reduction sequence, it is necessary to classify each file (e.g. dark, fflat, arcs, science) within every input observing night. This is done in py2dfdr by defining *observing_runs*, i.e., collections of nights that share a set of calibration files (namely, darks, bias and detector flats). 

The *archive.ArchiveObs* module provides a tool for classifying multiple observing nights into observing runs. The *input directory* must contain a set of nigth directories, using the naming convention: YYYYMMDD:
    |MySetOfNights
        |20220101
            ...FITS files...
        |20220102
            ...

All consecutive nights will be grouped into observing runs that extend up to two weeks of observations. The data withing each observing run will be organized according to the following scheme:

    |obs_run_0               (Observing run ID)
        |darks              (Dark files)
            |ccd1            
                |exptime1
                |exptime2
            |ccd2
                ...
        |bias               (Bias files)
            |ccd1
            |ccd2
        |fflats             (Detector flats) 
            |ccd1
                |grating1
                |grating2
                ...
            |ccd2
                ...
        |night_0            (Data corresponding to an individual night)
            |sci                 (cal. stars, skyflats, science...)
                |ccd1
                    |grating1
                    |grating2
                |ccd2
                    ...
            |arcs
                |ccd1
                    |grating1
                    |grating2
                |ccd2
            |fibreflat
                |ccd1
                    |grating1
                    |grating2
                |ccd2
        ...
        |night_n
    ...
    |obs_run_k

See examples in /tutorials for further instructions.

In order to prepare the data, use the following methods:
- Initialize the archive class
```commandline
from py2dfdr.archive import ArchiveObs

archive = ArchiveObs(path_to_data)
```
- Find the nights available and fetch the information for each file.
```commandline
archive.find_nights()
archive.get_night_info()
```
After this step, *archive.nigths* corresponds to a dictionary containing the information of each night.
- Finally, a new working folder containing the observing runs is created using the following method:
```commandline
archive.create_reduction_folder(output_path=new_working_dir_path, hard=False) 
```
If hard=False, it will create symbolic links instead of directly coping the original files.
A yaml file containing all the information for each observing run, 'obs_run_info.yml', will be saved.

## Observing run data reduction
Once the data is properly organised we can reduce each observing run.
```commandline
from reduce_obs_run import ReduceObsRun
### 
obsrunpath = '/media/pablo/toshiba-pab/reduce_koala_april/obs_run_0'
redOR = ReduceObsRun(
            obs_run_path=obsrunpath,
            ccds=['ccd_2'],
            dark_idx='koala_dark.idx',
            fibreflat_idx='koala_dark.idx',
            lflat_idx='koala_dark.idx',
            arcs_idx='koala_arcs.idx',
            object_idx='koala_reduce.idx')
```
For each calibration file you must provide an .idx configuration file.
Then, you can call the different routines to reduce sequentially each calibration/science .fits
```commandline
redOR.reduce_darks()
# redOR.get_master_darks()
redOR.reduce_lflats()
# redOR.get_master_lflats()
redOR.extract_tramlines()
# redOR.get_tlm_maps()
redOR.reduce_arcs()
# redOR.get_arcs()
redOR.reduce_fflats()
# redOR.get_fibreflats()
redOR.reduce_object()
```
During the process, some quality control plots will be generated.
If, for example, the observing run already contains the products after reducing all darks, it is possible to load them using 'get_master_dark'.

