#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 09:16:03 2022

@author: pablo
"""

from py2dfdr import KOALA_archive


# =============================================================================
# Create data reduction folder
# =============================================================================
archive = KOALA_archive('/media/pablo/toshiba-pab/KOALA_april')

archive.find_nights()

archive.create_observing_runs()

archive.get_night_info()

archive.create_reduction_folder(output_path='/media/pablo/toshiba-pab/reduce_koala_april')
