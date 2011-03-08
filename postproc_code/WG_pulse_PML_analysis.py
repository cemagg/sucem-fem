from __future__ import division
import pickle
from itertools import chain

import numpy as N
import scipy

from NewCode.Utilities import RMS, Struct

import analytical_WG_driver
from analytical_WG_driver import WG_test


def make_phase(tf, sl):
    return N.unwrap(N.angle(tf[sl]))

def merge_resses(resses):
    merged_res = {}
    for res in resses:
        for o, oresses in res.iteritems():
            print o
            if not o in merged_res: merged_res[o]=oresses
            else: merged_res[o].update(oresses)
    return merged_res


an_res = pickle.load(file('+WG_analytical_result_9s.pickle'))

ref_res_files = (
    #'+coupledb_brick_8_wg_1st-3rd_lumped.pickle',
    '+coupledb_brick_8_wg_1st-3rd_lumped_wglen_6.5.pickle',
    )

res_files = (
#     '+coupledb_brick_wg_PML_tmp.pickle',
#    '+brick_8_wg_PML_15c_1_m3.5.pickle',
    '+brick_8_wg_PML_tmp.pickle',
    )

# res_files = (
#     '+coupledb_brick_8_wg_1st-3rd_lumped.pickle',
#     )

ref_res = merge_resses(pickle.load(file(fn)) for fn in ref_res_files)
res = merge_resses(pickle.load(file(fn)) for fn in res_files)

f_min, f_max = 0.4375, 2.375
n_f_min, n_f_max = int(f_min/an_res.df), int(f_max/an_res.df)

freq_slice = slice(n_f_min, n_f_max)
freq_range = an_res.freqs[freq_slice]

an_res_phase = make_phase(an_res.transfer_fs, freq_slice)

num_analysis = {}

for tsr, order, div in chain(
    *[[(val,o,d) for d,val in r.items()] for o,r in res.items()]):
    ref_tsr = ref_res[order][div]
    WGN = WG_test()
    WGN.ts = tsr.ts_modeintg_n - ref_tsr.ts_modeintg_n
    WGN.calc_FFTs()
    refl_phase = make_phase(WGN.output_fs, freq_slice)
    refl_mag = N.abs(WGN.output_fs[freq_slice])
    refl_mag_RMS = RMS(refl_mag)
    dt = WGN.dt/div
    num_analysis[order, dt] = Struct(phase=refl_phase, refl_mag=refl_mag,
                                     refl_mag_RMS=refl_mag_RMS,
                                     orig_ts=ref_tsr.ts_modeintg_n,
                                     refl_ts=WGN.ts, dt=dt, an_dt=WGN.dt)
    
RMS_mag_errs = []
RMS_phase_errs = []
dts = []

f_c = WG_test.omega_c/2/N.pi
vg = lambda f: 2*N.pi*N.sqrt(1-(f_c/f)**2) if f > f_c else 0
Beta = lambda f: f*2*N.pi*N.sqrt(1-(f_c/f)**2) if f > f_c else 0

dts_o = dict((o, [dt for oo,dt in num_analysis.keys() if oo==o]) for o in res.keys())
max_dts = dict((o, max(dt for oo,dt in num_analysis.keys() if oo==o)) for o in res.keys())
min_dts = dict((o, min(dt for oo,dt in num_analysis.keys() if oo==o)) for o in res.keys())

print "(order, dt) : mag err"
for o in res.keys():
    mid_no = int(N.floor(len(dts_o[o])/2.))
    max_dt, mid_dt, min_dt = max_dts[o], sorted(dts_o[o])[mid_no], min_dts[o]
    print (o, max_dt), ' : ', N.log10(num_analysis[(o, max_dt)].refl_mag_RMS)
#     print (o, mid_dt), ' : ', N.log10(num_analysis[(o, mid_dt)].phase_err_RMS), \
#           ' ', N.log10(num_analysis[(o, mid_dt)].refl_mag_RMS)
    print (o, min_dt), ' : ', N.log10(num_analysis[(o, min_dt)].refl_mag_RMS)
