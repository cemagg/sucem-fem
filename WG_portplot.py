from __future__ import division
import numpy as N
import glob
import re
import os
import pickle

from NewCode.Utilities import Struct
from NewCode.PostProc.AnalyseTFSFPort import AnalyseTFSFPort
from NewCode.Analytical import WaveguidePhasor
from RunInfos.WGParms import WGParms
from NewCode.Consts import c0

line_types = ('k-.*', 'k--s', 'k--o', 'k-.D', 'k-^', 'k--x')
max_df = 1e6/c0                            # FFT freqency resolution
# res_dir = 'hcyl_pin_output'
#res_dir = 'multipost_output'
#res_fglob = 'hybmesh_wg_lech03fig7*.pickle'
res_dir = 'hcyl_pin_output'
res_fglob = '*_*wg_hcyl_pin_1*.pickle'
#res_fglob = '*brick*_wg_hcyl_pin_1*.pickle'
res_filenames = glob.glob(res_dir+os.sep+res_fglob)
info_regxp= re.compile(r'(.*)_\+?wg_.*_(h.*?)_o(\d*).*')
a=WGParms['WR90'].a
an_calc = WaveguidePhasor.TE01(a)
f_c = an_calc.k_c/2/N.pi
f_min, f_max = 1.05*f_c, 2*f_c
AnP = AnalyseTFSFPort(max_df)
AnP.set_freqrange(f_min, f_max)
# run_ffts = Struct()
# res_ffts = run_ffts
# all_dts = set()
# all_orders = set()

# for fn in res_filenames:
#     meshtype, meshsize, order = info_regxp.match(os.path.basename(fn)).groups()
#     meshname = meshtype+'_'+meshsize
#     order = int(order) 
#     try: mesh_ffts = run_ffts[meshname]
#     except KeyError: mesh_ffts = run_ffts[meshname] = dict()
#     try: order_ffts = mesh_ffts[order]
#     except KeyError: order_ffts = mesh_ffts[order] = dict()
#     res = pickle.load(file(fn))
#     AnP.set_dt(res.base_dt)
#     order_ffts[res.run_dt] = AnP.get_ScatParms(res.drv_ts, res.ts_modeintg1_n,
#                                        res.ts_modeintg2_n)
#     all_dts.add(res.run_dt)
#     all_orders.add(order)

base_dt = res.base_dt
all_dts = N.unique(all_dts)
dts_by_key = {}
for mshkey, mshval in res_ffts.iteritems():
    dts_by_key[mshkey] = {}
    for order, ordervals in mshval.iteritems():
        dts_by_key[mshkey][order] = N.unique(ordervals.keys())

all_orders = N.unique(all_orders)

import FEKO_WR90_hcyl_pin123 as FEK

ref_msh = 'hybmesh_h16'
ref_dt = 0.00012498321168252695
ref_o = 3
ref_linetype = 'k-'
ref_xlim = (9.9, 13.2)
ref_ylim = (-41, 1)
ref_xlabel = 'Frequency (GHz)'
ref_ylabel = '|S11| (dB)'
def do_o1_hybrid_plots():
    meshnames = ('hybmesh_h8',
                 'hybmesh_h16',
                 'hybmesh_h32',
                 'hybmesh_h64', )
    mesh_caps = ('h8', 'h16', 'h32', 'h64')
    pr = res_ffts[ref_msh][ref_o][ref_dt] 
    plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), ref_linetype, label='ref' )
    for i, (msh, cap) in enumerate(zip(meshnames, mesh_caps)):
        o = 1 ; dt = dts_by_key[msh][o][-1] ; pr = res_ffts[msh][o][dt]
        plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), line_types[i], label=cap, markevery=150 )
    xlim(ref_xlim)
    ylim(ref_ylim)
    xlabel(ref_xlabel)
    ylabel(ref_ylabel)
    legend(loc='lower left')
    
def do_o1_hexa_plots():
    meshnames = ( 'brick_bn_h8',
                  'brick_bn_h11',
                  'brick_bn_h16',
                  'brick_bn_h21',
                  'brick_bn_h32',
                  #'brick_bn_h48',
                  'brick_bn_h64',
                  'brick_bn_h95',
                  'brick_bn_h128')
    mesh_caps = ('h8', 'h11', 'h16', 'h21', 'h32', #'h48',
                 'h64', 'h95', 'h128')
    pr = res_ffts[ref_msh][ref_o][ref_dt] 
    plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), ref_linetype, label='ref' )
    for i, (msh, cap) in enumerate(zip(meshnames, mesh_caps)[3:]):
        o = 1 ; dt = dts_by_key[msh][o][-1] ; pr = res_ffts[msh][o][dt]
        plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), line_types[i], label=cap, markevery=150 )
    xlim(ref_xlim)
    ylim(ref_ylim)
    xlabel(ref_xlabel)
    ylabel(ref_ylabel)
    legend(loc='lower left')

def do_o2_hybrid_plots():
    meshnames = ('hybmesh_h4',
                 'hybmesh_h8',
                 'hybmesh_h16',
                 'hybmesh_h32',
                  )
    mesh_caps = ('h4', 'h8', 'h16', 'h32')
    o = 2 
    pr = res_ffts[ref_msh][ref_o][ref_dt] 
    plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), ref_linetype, label='ref' )
    for i, (msh, cap) in enumerate(zip(meshnames, mesh_caps)[0:-1]):
        pr = res_ffts[msh][o][ref_dt]
        plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), line_types[i], label=cap, markevery=150 )
    xlim(ref_xlim)
    ylim(ref_ylim)
    xlabel(ref_xlabel)
    ylabel(ref_ylabel)
    legend(loc='lower left')


def do_o2_hybrid_refined_plots():
    meshnames = ('hybmesh_h3_refined1',
                 'hybmesh_h4_refined1m',
                 'hybmesh_h8_refined1',
                 'hybmesh_h16',
                 'hybmesh_h32',
                  )
    mesh_caps = ('h3_refined', 'h4_refined', 'h8_refined', 'h16', 'h32')
    o = 2 
    pr = res_ffts[ref_msh][ref_o][ref_dt] 
    plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), ref_linetype, label='ref' )
    for i, (msh, cap) in enumerate(zip(meshnames, mesh_caps)[1:-2]):
        pr = res_ffts[msh][o][ref_dt]
        plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), line_types[i], label=cap, markevery=150 )
    xlim(ref_xlim)
    ylim(ref_ylim)
    xlabel(ref_xlabel)
    ylabel(ref_ylabel)
    legend(loc='lower left')

def do_o3_hybrid_plots():
    meshnames = ('hybmesh_h3',
                 'hybmesh_h4',
                 'hybmesh_h8',
                  )
    mesh_caps = ('h3', 'h4', 'h8',)
    o = 3 
    pr = res_ffts[ref_msh][ref_o][ref_dt] 
    plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), ref_linetype, label='ref' )
    for i, (msh, cap) in enumerate(zip(meshnames, mesh_caps)):
        pr = res_ffts[msh][o][ref_dt]
        plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), line_types[i], label=cap, markevery=150 )
    xlim(ref_xlim)
    ylim(ref_ylim)
    xlabel(ref_xlabel)
    ylabel(ref_ylabel)
    legend(loc='lower left')


def do_o3_hybrid_refined_plots():
    meshnames = ('hybmesh_h3_refined1',
                 'hybmesh_h4_refined1m',
                 'hybmesh_h8_refined1',
                  )
    mesh_caps = ('h3_refined', 'h4_refined', 'h8_refined')
    o = 3
    pr = res_ffts[ref_msh][ref_o][ref_dt] 
    plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), ref_linetype, label='ref' )
    for i, (msh, cap) in enumerate(zip(meshnames, mesh_caps)):
        pr = res_ffts[msh][o][ref_dt]
        plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), line_types[i], label=cap, markevery=150 )
    xlim(ref_xlim)
    ylim(ref_ylim)
    xlabel(ref_xlabel)
    ylabel(ref_ylabel)
    legend(loc='lower left')

def do_o3_o2_hybrid_refined():
    meshnames = (
        'redhybmesh_h4_refined1m',
        'redhybmesh_h4_refh8',
        'hybmesh_h8',
                  )
    mesh_caps = ('o3_h4_refined_o2', 'o3_h4_refh8_o2', 'o2_h8')
    orders = [3, 3, 2]
    pr = res_ffts[ref_msh][ref_o][ref_dt] 
    plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), ref_linetype, label='ref' )
    for i, (msh, cap, o)  in enumerate(zip(meshnames, mesh_caps, orders)[1:]):
        pr = res_ffts[msh][o][ref_dt]
        plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), line_types[i], label=cap, markevery=150 )
    xlim(ref_xlim)
    ylim(ref_ylim)
    xlabel(ref_xlabel)
    ylabel(ref_ylabel)
    legend(loc='lower left')


"""
plot(FEK.fr_ref/1e9, FEK.pin1_ref_S11, line_types[0], label='FEKO ref')
msh = 'hybmesh_h16' ; o = 3 ; dt = all_dts[2] ; pr = res_ffts[msh][o][dt] ; plot(pr.fr*c0/1e9, 20*N.log10(pr.S11), line_types[1], label='hybrid ref' )

xlim(8, 13.5)
ylim(-40, 0)
"""
"""
msh = 'hybmesh_h16' ; o = 3 ; dt = all_dts[3] ; pr = res_ffts[msh][o][dt] ; plot(pr.fr*c0, 20*N.log10(pr.S11), label='m:%s, o:%d, dt:%f' % ( msh, o, dt))
"""
