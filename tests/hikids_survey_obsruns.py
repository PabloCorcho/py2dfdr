from py2dfdr.archive import ArchiveObs

# Path to folder containing all observing nights (e.g. 220112...)
path_to_data = '/Users/pcorcho/obs_data/KOALA/raw'
archive = ArchiveObs(path_to_data)
archive.find_nights()
archive.get_night_info()
archive.create_ori()
archive.create_reduction_folder(
    output_path='/Users/pcorcho/obsruns', hard=False)