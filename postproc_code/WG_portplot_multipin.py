from __future__ import division
import numpy as N
import glob
import re
import os
import pickle

import scipy.signal
from NewCode.Utilities import Struct, partial
from NewCode.PostProc.AnalyseTFSFPort import AnalyseTFSFPort
from NewCode.Analytical import WaveguidePhasor
from RunInfos.WGParms import WGParms
from NewCode.Consts import c0
from NewCode.PostProc import AnalyseTFSFPort
reload(AnalyseTFSFPort)
AnalyseTFSFPort = AnalyseTFSFPort.AnalyseTFSFPort

def tukey(n, a=0.25):
    p0 = int(N.ceil(a/2*(n-1)))
    p1 = int(n - N.floor(a/2*(n-1)))
    k0 = N.arange(0, p0)
    k2 = N.arange(p1, n)
    return N.hstack([1/2*(1 + N.cos(2*N.pi/a*k0/(n-1) - N.pi)),
                     N.ones(p1-p0, dtype=N.float64),
                     1/2*(1 + N.cos(2*N.pi/a*(1-k2/(n-1)) - N.pi))])

line_types = ('k-o', 'k-.s', 'k--D', 'k--o', 'k-.D', 'k-^', 'k--x')
max_df = 1e6/c0                            # FFT freqency resolution
# res_dir = 'hcyl_pin_output'
# res_fglob = '*_pin_1*.pickle'
res_dir = 'multipost_output'
fig = 'fig4'
#fig = 'fig7'
#res_fglob = '*hybmesh_wg_lech03fig4r4.41*o[23]*.pickle'
res_fglob = '*hybmesh_wg_lech03%s*o[23]*.pickle' % fig
res_filenames = glob.glob(res_dir+os.sep+res_fglob)
info_regxp= re.compile(r'(.*)_\+?wg_.*_(h.*?)_o(\d*).*')
runfac_regxp = re.compile(r'.*_runfac(\d*\.?\d*).*')
a=WGParms['WR90'].a
an_calc = WaveguidePhasor.TE01(a)
f_c = an_calc.k_c/2/N.pi
#f_min, f_max = 1.05*f_c, 2*f_c
#f_min, f_max = N.array([10e9, 10.4e9])/c0
f_minmax = Struct(fig7=[1.05*f_c, 1.99*f_c],
                  fig4=N.array([9e9, 11.e9])/c0)
                               
f_min, f_max = f_minmax[fig]
AnP = AnalyseTFSFPort(max_df)
#AnP_w = AnalyseTFSFPort(max_df, partial(tukey, a=0.4))
#AnP_w = AnalyseTFSFPort(max_df, scipy.signal.hann)
AnP_w = AnalyseTFSFPort(max_df, scipy.signal.hamming)
#AnP_w = AnalyseTFSFPort(max_df, scipy.signal.nuttall)
AnP.set_freqrange(f_min, f_max)
AnP_w.set_freqrange(f_min, f_max)
res_ffts = Struct()

all_dts = set()
all_orders = set()

for fn in res_filenames:
    meshtype, meshsize, order = info_regxp.match(os.path.basename(fn)).groups()
    meshname = meshtype+'_'+meshsize
    runfac_m = runfac_regxp.match(fn)
    if runfac_m: meshname += '_rf'+runfac_m.groups()[0]
    order = int(order) 
    try: mesh_ffts = res_ffts[meshname]
    except KeyError: mesh_ffts = res_ffts[meshname] = dict()
    try: order_ffts = mesh_ffts[order]
    except KeyError: order_ffts = mesh_ffts[order] = dict()
    res = pickle.load(file(fn))
    AnP.set_dt(res.base_dt)
    AnP_w.set_dt(res.base_dt)
    order_ffts[res.run_dt] = Struct(nw=AnP.get_ScatParms(res.drv_ts, res.ts_modeintg1_n,
                                                         res.ts_modeintg2_n),
                                    w=AnP_w.get_ScatParms(res.drv_ts, res.ts_modeintg1_n,
                                                          res.ts_modeintg2_n))
    all_dts.add(res.run_dt)
    all_orders.add(order)

all_dts = N.unique(all_dts)
all_orders = N.unique(all_orders)
dts_by_key = {}
for mshkey, mshval in res_ffts.iteritems():
    dts_by_key[mshkey] = {}
    for order, ordervals in mshval.iteritems():
        dts_by_key[mshkey][order] = N.unique(ordervals.keys())

import FEKO_WR90_lech03_fig4 as FEK4
import FEKO_WR90_lech03_fig7 as FEK7
ref_dt = 0.00012498321168252695
ref_dt_fig4 = 0.00024996642336505391
ref_xlim_fig4 = (10.1, 10.4)
ref_ylim_fig4 = (-40, 0)
ref_xlim_fig7 = (8, 13.2)
ref_ylim_fig7 = (-25, 0)
from pylab import xlim, ylim, legend, plot

def do_plots_fig4(sl=slice(10000), win='nw'):
    meshnames = (
        'hybmesh_h4_refined1_rf4.000000',
        'hybmesh_h4_refined1_rf8.000000',
        'hybmesh_h8_refined1_rf8.000000',
        'hybmesh_h8_refined2_rf8.000000',
        'hybmesh_h8_refined2_rf8.000000',
        'redhybmesh_h8_refined2_rf4.000000')
    mesh_caps = (
        'o3_h4_r1_rf4',#0
        'o3_h4_r1_rf8',#1
        'o3_h8_r1_rf8',#2
        'o2_h8_r2_rf8',#3
        'o3_h8_r2_rf8',#4
        'o3_h8_r2_rf4_o2')
    orders = [3,3,3,2,3,3]
    plot(FEK4.fr_ref2/1e9, FEK4.r441_S11_ref2, '--', label='FEKO')
    for msh, cap, o in zip(meshnames, mesh_caps, orders)[sl]:
        pr = res_ffts[msh][o][ref_dt][win]
        plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), '--', label='hybrid') #label=cap)
    xlim(ref_xlim_fig4)
    ylim(ref_ylim_fig4 )
    
def do_plots_fig7(sl=slice(10000), win='nw'):
    meshnames = (
        'hybmesh_h3_refined1',
        'hybmesh_h4_refined1',
        'hybmesh_h4_refined1',
        'hybmesh_h4_refh8',
        'hybmesh_h4_refh8r1',
        'hybmesh_h8',
        'redhybmesh_h4_refh8',
        'redhybmesh_h4_refh8r1')
    mesh_caps = (
        'o3_h3_r1',
        'o2_h4_r1',
        'o3_h4_r1',
        'o3_h4_rh8',
        'o3_h4_rh8r1',
        'o3_h8',
        'o3_h4_rh8_o2',
        'o3_h4_rh8r1_o2')
    orders = [3,2,3,3,3,3,3,3]
    plot(FEK7.fr/1e9, FEK7.S11, label='fek')
    for msh, cap, o in zip(meshnames, mesh_caps, orders)[sl]:
        pr = res_ffts[msh][o][ref_dt][win]
        plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), label=cap)
        xlim(ref_xlim_fig7)
        ylim(ref_ylim_fig7)

# f_min, f_max = 0.875*f_c, 2.2*f_c
# res_ts = resses[order][dt_div].ts_modeintg2_n
# res_ts1 = resses[order][dt_div].ts_modeintg1_n
# n_ts = len(res_ts)
# #n = 2**int(N.ceil(N.log2(n_ts)))
# n = 2**16
# assert(n > n_ts)
# f_nyq = 1/2/dt/dt_div
# df = f_nyq/n*2
# res_fs = scipy.fft(res_ts, n=n)
# res_fs1 = scipy.fft(res_ts1, n=n)
# drv_fs = scipy.fft(drv_ts, n=n)
# n_f_min, n_f_max = int(N.floor(f_min/df)), int(N.ceil(f_max/df))
# tf = (res_fs/drv_fs)[n_f_min:n_f_max]
# tf1 = ((res_fs1-drv_fs)/drv_fs)[n_f_min:n_f_max]
# tf1_tf = (res_fs1/drv_fs)[n_f_min:n_f_max]
# fr = N.arange(n_f_min, n_f_max)*df*Consts.c0
# tr = N.arange(0, no_steps+1, dt_div)*dt


