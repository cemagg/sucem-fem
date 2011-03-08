"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
from itertools import izip
import os, sys
import pickle
import random
import numpy as N
import scipy
#
# Local Imports
#
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import Waveforms
import NewCode.eMAGUSImport as eMAGUSImport
from NewCode.Runners.HybridWaveGuide2Port import HybridWaveGuide2PortRect

from analytical_WG_driver import WGT

h = 1/4.
a,b = 1., 0.25
x_c = 1/3
z_c = 5/8
d=1/5
tet_bot_left = N.array([x_c-d/2, 0, z_c-d/2], N.float64)
tet_top_right = N.array([x_c+d/2, b, z_c+d/2], N.float64)


dt_div = 4
dt = 1/16/dt_div
order = 2
l_sc = h
z_inc = 0
z_measure = WGT.test_z + z_inc
drv_fun = WGT.discrete_drv_fn
#runtime = WGT.t_final
runtime = 125
source_runtime = None

z_PML_goal = z_measure
z_PML1_goal = -1/2
assert (z_PML1_goal <= -l_sc)
no_PML_cells = 30
source_no_PML_cells = 30

eMAGUSImport.init('workspace')
tet_mesh = eMAGUSImport.get_mesh()

print 'order: ', order
TRS = HybridWaveGuide2PortRect(a,b,h,order,drv_fun, no_PML_cells)
TRS.log_divisor = dt_div
TRS.set_dims(z_inc, z_inc, z_measure)
TRS.init_mesh()
TRS.init_hyb_geom(tet_bot_left, tet_top_right)
TRS.init_tetmesh(tet_mesh)
TRS.init_systems()
TRS.setupSource(source_runtime=source_runtime, source_no_PML_cells=source_no_PML_cells)
TRS.set_dt(dt)
no_steps = int(N.ceil(runtime/dt))
TRS.runSteps(no_steps)
resses = {order:{dt_div:TRS.getResult()}}
drv_ts = N.array(TRS.inc_modeintg[::dt_div], N.float64)

f_min, f_max = 0.4375, 2.375
res_ts = resses[order][dt_div].ts_modeintg2_n
res_ts1 = resses[order][dt_div].ts_modeintg1_n
n_ts = len(res_ts)
#n = 2**int(N.ceil(N.log2(n_ts)))s
n = 2**14
assert(n > n_ts)
df = WGT.f_nyq/n*2
res_fs = scipy.fft(res_ts, n=n)
res_fs1 = scipy.fft(res_ts1, n=n)
drv_fs = scipy.fft(drv_ts, n=n)
n_f_min, n_f_max = int(N.floor(f_min/df)), int(N.ceil(f_max/df))
tf = (res_fs/drv_fs)[n_f_min:n_f_max]
tf1 = ((res_fs1-drv_fs)/drv_fs)[n_f_min:n_f_max]
tf1_tf = (res_fs1/drv_fs)[n_f_min:n_f_max]
fr = N.arange(n_f_min, n_f_max)*df
tr = N.arange(0, no_steps+1, dt_div)*dt

from NewCode.Analytical import WaveguidePhasor

an_calc = WaveguidePhasor.TE01(a)

z_tf = z_measure - z_inc 
an_tf = N.array([an_calc.TF(f*2*N.pi, z_tf) for f in fr], N.complex128)
