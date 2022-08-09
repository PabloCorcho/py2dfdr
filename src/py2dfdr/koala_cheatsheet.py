"""
This script contains the optimal calibration parameters for matching the correct reduction
files.

DISCLAIMER:
 Exposure times and correct arc lamps have only been tested with the low resolution gratings. The rest of the values
 might be subjected to errors and communication with the instrument scientist is highly recommended before reducing the
 data.
"""

# =============================================================================
# Recommended fibre flat exposure times for each grating
# =============================================================================

fflat_time = {
    'ccd_1':
        {
         '1700B': 45,
         '3200B': 45,
         '580V': 45.,
         '1500V': 30,
         '2500V': 30},
    'ccd_2':
        {
         '385R': 4.,
         '1000R': 4.,
         '2000R': 4.,
         '1000I': 4.,
         '1700I': 10,
         '1700D': 10,
         }
                }

# =============================================================================
# Recommended long-slit (detector) flats exposure times for each grating
# =============================================================================

lflat_time = {
    'ccd_1':
        {
         '1700B': 45,
         '3200B': 45,
         '580V': 45.,
         '1500V': 30,
         '2500V': 30},
    'ccd_2':
        {
         '385R': 2.,
         '1000R': 45.,
         '2000R': 4.,
         '1000I': 4.,
         '1700I': 10,
         '1700D': 10,
         }
                }

# =============================================================================
# Recommended arcs and exposure times for each grating
# =============================================================================

arc_time = {
    'ccd_1':
        {
         '1700B': ['CuAr', 120.],
         '3200B': ['CuAr', 120.],
         '580V': ['CuAr', 40.],
         '1500V': ['CuAr', 120.],
         '2500V': ['CuAr', 120.]
         },
    'ccd_2':
        {
         '385R': ['CuAr', 0.1],
         '1000R': ['CuAr', 3.0],
         '2000R': ['CuAr', 10.0],
         '1000I': ['CuAr', 3.0],
         '1700I': ['CuAr', 10.0],
         '1700D': ['CuAr', 10.0],
        }
        }