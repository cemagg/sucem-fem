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


an_res = pickle.load(file('+WG_analytical_result_9s.pickle'))
# res_files = (
#     '+newmark_hybrid_1st_8_tmp.pickle',
#     '+newmark_hybrid_2nd_8_tmp.pickle',
#     '+newmark_hybrid_3rd_8.pickle',
#     '+newmark_hybrid_4th_8.pickle',
#    '+newmark_wg_hybmesh_tmp.pickle',
#    '+newmark_wg_hybmesh-2.pickle',
#    '+newmark_wg_hybmesh_tmp-4.pickle',
#    '+newmark_wg_hybmesh_tmp-8.pickle',
#    '+newmark_wg_hybmesh_tmp-16.pickle',
#    '+newmark_hybrid_wg_hybmesh_tmp.pickle',
#    '+newmark_hybrid_wg_hybmesh-2.pickle',
#    '+newmark_hybrid_wg_hybmesh-4.pickle',
#    '+newmark_hybrid_wg_hybmesh-8.pickle',
#    '+newmark_hybrid_wg_hybmesh-16.pickle',
#              )
# res_files = (
#     '+newmark_leapfrog_hybrid_1st_8.pickle',
#     '+newmark_leapfrog_hybrid_2nd_8.pickle',
#     '+newmark_leapfrog_hybrid_3rd_8_tmp.pickle',
#     '+newmark_leapfrog_hybrid_4th_8_tmp.pickle',
#              )

# res_files = (
#     '+coupledb_brick_wg_1st-4th_consistent.pickle',
#     )

# res_files = (
#     '+newmark_hybrid_wg_hybmesh-tmp.pickle',
#     )

# res_files = (
#     '+newmark_leapfrog_hybrid_tmp.pickle', 
#     )

# res_files = (
#     '+newmark_leapfrog_hybrid_hybmesh-tmp.pickle', 
#     )

# res_files = (
#     '+brick_4_wg_PML_tmp.pickle',
#     )

res_files = (
    '+brick_4_wg_tfsf_tmp.pickle',
    )
# res_files = (
#     '+newmark_pml_hybrid_tmp.pickle', 
#     )


# res_files = (
#     '+newmark_hybrid_wg_hybmesh_start-tmp.pickle',
#     )

# res_files = (
#     '+newmark_pml_wg_hybmesh-tmp.pickle',
#     )

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
    if hasattr(tsr, 'drv_ts'):
        WGN.calc_FFTs_with_drv_ts(tsr.drv_ts)
        print "own drv_ts!!!!!!"
    else: WGN.calc_FFTs()

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

#pickle.dump(num_analysis, file('+WG_pulse_analysis_newmark_8.pickle', 'w'))

dts_o = dict((o, [dt for oo,dt in num_analysis.keys() if oo==o]) for o in res.keys())
max_dts = dict((o, max(dt for oo,dt in num_analysis.keys() if oo==o)) for o in res.keys())
min_dts = dict((o, min(dt for oo,dt in num_analysis.keys() if oo==o)) for o in res.keys())

print "(order, dt) : phase err   mag err"
for o in res.keys():
    mid_no = int(N.floor(len(dts_o[o])/2.))
    max_dt, mid_dt, min_dt = max_dts[o], sorted(dts_o[o])[mid_no], min_dts[o]
    print (o, max_dt), ' : ', N.log10(num_analysis[(o, max_dt)].phase_err_RMS), \
          ' ', N.log10(num_analysis[(o, max_dt)].mag_err_RMS)
#     print (o, mid_dt), ' : ', N.log10(num_analysis[(o, mid_dt)].phase_err_RMS), \
#           ' ', N.log10(num_analysis[(o, mid_dt)].mag_err_RMS)
    print (o, min_dt), ' : ', N.log10(num_analysis[(o, min_dt)].phase_err_RMS), \
          ' ', N.log10(num_analysis[(o, min_dt)].mag_err_RMS)
