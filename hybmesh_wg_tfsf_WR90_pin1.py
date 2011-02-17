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
from NewCode.Runners.HybridWaveGuide2Port import HybridWaveGuide2PortPECRect
from NewCode.Analytical import WaveguidePhasor
from NewCode import Consts

a,b = 22.86e-3, 10.16e-3
#a = 22.86e-3
#b=2/5*a
an_calc = WaveguidePhasor.TE01(a)

f_c = an_calc.k_c/2/N.pi
lam_c = 1/f_c
h = (lam_c/2)/5.
order = 2
dt_div = 6
# dt_div = 1 corresponds to just-above stability dt for a/5 1st order element
base_dt = lam_c/2/5*Consts.lumped_stability_factors[1]/1.1
dt_div_min = int(N.ceil(base_dt/(h*Consts.lumped_stability_factors[order])))
assert(dt_div >= dt_div_min)            # Not conservative enough since cells can be
                                        # smaller in y direction
dt = base_dt/dt_div 

pin_z_cent = a*(1/2 + 1/5)
pin_x_cent = a/5
pin_d = a/5

l_sc = h
z_inc = 0
z_measure = (lam_c/2)*7/5
drv_fun = Waveforms.get_gausspulse(fc=1.5*f_c, bw=1.1,tpr=-60)

l_sc = h
z_inc = 0

runtime = 62.5/f_c
source_runtime = runtime

no_PML_cells = 15
source_no_PML_cells = no_PML_cells


wiggle_room = N.array([h/2, 0, h/2], N.float64)

tet_bot_left = N.array([a/2-pin_d/2+pin_x_cent, 0, pin_z_cent-pin_d/2],
                   N.float64) - wiggle_room
tet_top_right = N.array([a/2+pin_d/2+pin_x_cent, b, pin_z_cent+pin_d/2],
                    N.float64) + wiggle_room

eMAGUSImport.init('workspace')
tet_mesh = eMAGUSImport.get_mesh()

print 'order: ', order
TRS = HybridWaveGuide2PortPECRect(a,b,h,order,drv_fun, no_PML_cells)
TRS.log_divisor = dt_div
TRS.set_dims(z_inc, z_inc, z_measure)
TRS.init_mesh()
TRS.init_hyb_geom(tet_bot_left, tet_top_right)
TRS.init_tetmesh(tet_mesh)
TRS.init_hybrid_PEC(PEC_label=1001)
TRS.init_systems()
TRS.setupSource(source_runtime=source_runtime, source_no_PML_cells=source_no_PML_cells)
TRS.set_dt(dt)
no_steps = int(N.ceil(runtime/dt))
TRS.runSteps(no_steps)
resses = {order:{dt_div:TRS.getResult()}}
drv_ts = N.array(TRS.inc_modeintg[::dt_div], N.float64)

f_min, f_max = 0.875*f_c, 2.2*f_c
res_ts = resses[order][dt_div].ts_modeintg2_n
res_ts1 = resses[order][dt_div].ts_modeintg1_n
n_ts = len(res_ts)
#n = 2**int(N.ceil(N.log2(n_ts)))
n = 2**14
assert(n > n_ts)
f_nyq = 1/2/dt/dt_div
df = f_nyq/n*2
res_fs = scipy.fft(res_ts, n=n)
res_fs1 = scipy.fft(res_ts1, n=n)
drv_fs = scipy.fft(drv_ts, n=n)
n_f_min, n_f_max = int(N.floor(f_min/df)), int(N.ceil(f_max/df))
tf = (res_fs/drv_fs)[n_f_min:n_f_max]
tf1 = ((res_fs1-drv_fs)/drv_fs)[n_f_min:n_f_max]
tf1_tf = (res_fs1/drv_fs)[n_f_min:n_f_max]
fr = N.arange(n_f_min, n_f_max)*df
tr = N.arange(0, no_steps+1, dt_div)*dt

# z_tf = z_measure - z_inc 
# an_tf = N.array([an_calc.TF(f*2*N.pi, z_tf) for f in fr], N.complex128)
