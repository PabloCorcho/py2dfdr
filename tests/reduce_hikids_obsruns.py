from py2dfdr import reduce_obs_run as red_obs

import time
tstart = time.time()
obsrunpath = '/Users/pcorcho/obs_data/KOALA/obsruns/obs_run_0'
redOR = red_obs.ReduceObsRun(
        obs_run_path=obsrunpath,
        ccds=[
              'ccd_1',
              'ccd_2'
              ],
        dark_idx='koala_dark.idx',
        lflat_idx='koala_dark.idx',
        fibreflat_idx='koala_fflat.idx',
        arcs_idx='koala_arcs.idx',
        object_idx='koala_reduce.idx')

# redOR.reduce_darks(timeout=300)
# redOR.get_master_darks()

# redOR.reduce_lflats(timeout=300)
# redOR.get_master_lflats()

redOR.extract_tramlines(timeout=300)
# redOR.get_master_tlm()

redOR.reduce_arcs(timeout=600)
# redOR.get_arcs()

# TODO: This should only reduce the master fibre flats selected during tramline extraction
redOR.reduce_fflats(timeout=600)
redOR.get_fibreflats()

# redOR.reduce_object(timeout=900)







tend = time.time()
print('\n\n ### Elapsed time (hrs): ', (tend - tstart) / 3600)

# Mr Krtxo \(ﾟ▽ﾟ)/
