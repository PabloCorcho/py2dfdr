#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 08:27:54 2022

@author: pablo
"""

"""
This script contains the optical parameters for matching the right reduction
files
"""

# =============================================================================
# Recommended exposure times
# =============================================================================

fflat_time = {
    'ccd_1':
        {'580V': 45.},
    'ccd_2':
        {'1000R': 45.,
         '385R': 4.}
                }

lflat_time = {
    'ccd_1':
        {'580V': 45.},
    'ccd_2':
        {'1000R': 45.,
         '385R': 2.}
                }

arc_time = {
    'ccd_1':
        {'580V': ['CuAr', 40.]},
    'ccd_2':
        {'385R': ['CuAr', 0.1],
         '1000R': ['CuAr', 3.0]}
        }