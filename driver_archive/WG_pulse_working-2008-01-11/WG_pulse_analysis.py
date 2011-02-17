from __future__ import division
import pickle
from itertools import chain

import numpy as N
import scipy
from scipy.interpolate import UnivariateSpline

from NewCode.Utilities import RMS, Struct

import analytical_WG_driver
from analytical_WG_driver import WG_test

interp_tol = 1e-7

def make_phase_spline(tf, sl):
    return UnivariateSpline(an_res.freqs[sl],
                            N.unwrap(N.angle(tf[sl])),
                            s=interp_tol)
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

# an_res = pickle.load(file('+WG_analytical_result_9s_wrong_drv.pickle'))
# res_files = ('+newmark_wg_1st_div1-26.pickle',
#              '+newmark_wg_2nd_div1-26.pickle',
#              '+newmark_wg_3rd_div1-26.pickle',
#              '+newmark_wg_4th_div1-26.pickle',
#              )

an_res = pickle.load(file('+WG_analytical_result_9s_new.pickle'))
# res_files = ('+newmark_wg_direch_1st_div1-128.pickle',
#              '+newmark_wg_direch_2nd_div1-128.pickle',
#              '+newmark_wg_direch_3rd_div1-128.pickle',
#              '+newmark_wg_direch_4th_div1-128.pickle',
#     )

# res_files = ('+coupledb_wg_direch_1st_div1-128.pickle',
#              '+coupledb_wg_direch_2nd_div1-128.pickle',
#              '+coupledb_wg_direch_3rd_div1-128.pickle',
#              '+coupledb_wg_direch_4th_div1-128.pickle',
#              )
res_files = ('+coupleda_wg_direch_1st_div1-128.pickle',
             '+coupleda_wg_direch_2nd_div1-128.pickle',
             '+coupleda_wg_direch_3rd_div1-128.pickle',
             '+coupleda_wg_direch_4th_div1-128.pickle',
             )
# res_files = ('+coupledb_wg_direch_3rd_9_div1-256.pickle',
#              '+coupledb_wg_direch_4th_9_div1-256.pickle',
#              )
# res_files = ('+coupledb_wg_direch_2nd_9_div1-256_newBCs.pickle',
#              '+coupledb_wg_direch_3rd_9_div1-256_newBCs.pickle',
#              '+coupledb_wg_direch_4th_9_div1-256_newBCs.pickle',
#              )

# res_files = ('+coupledb_wg_direch_3rd_11_div1-256_newBCs.pickle',
#              '+coupledb_wg_direch_2nd_11_div1-256_newBCs.pickle',
#              )
# res_files = ('+coupledb_wg_direch_1st_15_div1-256.pickle',
#              '+coupledb_wg_direch_2nd_15_div1-256.pickle',
#              )
             
# res_files = ('+coupledb_wg_direch_1st_15_div1-256.pickle',
#              '+coupledb_wg_direch_2nd_15_div1-256_newBCs.pickle',
#              )


res = merge_resses(pickle.load(file(fn)) for fn in res_files)

f_min, f_max = 0.4375, 2.375
extra_spline_pts = 3
n_f_min, n_f_max = int(f_min/an_res.df), int(f_max/an_res.df)

freq_slice = slice(n_f_min, n_f_max)
freq_range = an_res.freqs[freq_slice]

spline_slice = slice(n_f_min-extra_spline_pts,n_f_max+extra_spline_pts)
an_res_phase_spline = make_phase_spline(an_res.transfer_fs, spline_slice)
an_res_phase = make_phase(an_res.transfer_fs, freq_slice)

num_analysis = {}

for tsr, order, div in chain(
    *[[(val,o,d) for d,val in r.items()] for o,r in res.items()]):
    WGN = WG_test()
    WGN.ts = tsr.ts_modeintg_n
    WGN.calc_FFTs()
    num_phase = make_phase(WGN.transfer_fs, freq_slice)
    num_phase_spline = make_phase_spline(WGN.transfer_fs, spline_slice)
    phase_err = N.abs(num_phase-an_res_phase)
    phase_err_RMS = RMS(phase_err)
    mag_err = N.abs(N.abs(WGN.transfer_fs[freq_slice]) 
                    - N.abs(an_res.transfer_fs[freq_slice]))
    mag_err_RMS = RMS(mag_err)
    dt = WGN.dt/div
    num_analysis[order, dt] = Struct(phase=num_phase, phase_err=phase_err,
                                     phase_err_RMS=phase_err_RMS,
                                     phase_spline=num_phase_spline,
                                     transfer_fs=WGN.transfer_fs,
                                     mag_err=mag_err, mag_err_RMS=mag_err_RMS,
                                     dt=dt)
    
RMS_mag_errs = []
RMS_phase_errs = []
dts = []
for o,d in sorted(num_analysis.keys()):
    na = num_analysis[o,d]
    #    plot(freq_range, N.log10(na.mag_err), label=str(d))
    RMS_mag_errs.append(na.mag_err_RMS)
    RMS_phase_errs.append(na.phase_err_RMS)
    dts.append(na.dt)
#legend(loc=0)

# figure()
# plot(N.log10(dts), N.log10(RMS_mag_errs))
# figure()
# plot(N.log10(dts), N.log10(RMS_phase_errs))

f_c = WG_test.omega_c/2/N.pi
vg = lambda f: 2*N.pi*N.sqrt(1-(f_c/f)**2) if f > f_c else 0
Beta = lambda f: f*2*N.pi*N.sqrt(1-(f_c/f)**2) if f > f_c else 0

pickle.dump(num_analysis, file('+WG_pulse_analysis_coupleda.pickle', 'w'))
