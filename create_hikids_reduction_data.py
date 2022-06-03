#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from py2dfdr import KOALA_archive

path_to_data = '/media/pablo/toshiba-pab/KOALA_april'  # Path to folder containing all observing nights (e.g. 220112...)
archive = KOALA_archive(path_to_data)

archive.find_nights()
archive.create_observing_runs()
archive.get_night_info()
archive.create_reduction_folder(output_path=path_to_data.replace('KOALA_april', 'reduce_koala_april'))