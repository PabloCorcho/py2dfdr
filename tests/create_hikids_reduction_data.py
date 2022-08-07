
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from py2dfdr import KOALA_archive

# Path to folder containing all observing nights (e.g. 220112...)
path_to_data = '/home/pablo/Research/obs_data/HI-KIDS/raw/mar2022'
archive = KOALA_archive(path_to_data)
archive.find_nights()
archive.get_night_info()
archive.create_observing_runs()
archive.create_reduction_folder(
    output_path=path_to_data.replace('mar2022', 'mar2022_obsruns'), hard=False)