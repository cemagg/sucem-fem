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
from NewCode.Runners.WaveGuide2Port import WaveGuide2Port
from NewCode.Analytical import WaveguidePhasor
from NewCode import Consts

a,b = 22.86e-3, 10.16e-3
#a = 22.86e-3
#b=2/5*a
an_calc = WaveguidePhasor.TE01(a)

f_c = an_calc.k_c/2/N.pi
lam_c = 1/f_c
h = a/4.
order = 2
dt_div = 32
base_dt = a/3*Consts.lumped_stability_factors[1]/1.1
#base_dt = lam_c/2/5*Consts.lumped_stability_factors[1]/1.1
dt_div_min = int(N.ceil(base_dt/(h*Consts.lumped_stability_factors[order])))
assert(dt_div >= dt_div_min)
dt = base_dt/dt_div 
l_sc = h
z_inc = 0
z_measure = lam_c/2
drv_fun = Waveforms.get_gausspulse(fc=1.5*f_c, bw=1.1,tpr=-60)

runtime = 62.5/f_c
source_runtime = runtime

z_PML_goal = z_measure
z_PML1_goal = -1/2
assert (z_PML1_goal <= -l_sc)
no_PML_cells = 25
source_no_PML_cells = no_PML_cells


print 'order: ', order
TRS = WaveGuide2Port(a,b,h,order,drv_fun, no_PML_cells)
TRS.log_divisor = dt_div
#TRS.set_dims(z_inc, z_PML1_goal, z_measure)
TRS.set_dims(z_inc, z_inc, z_measure)
TRS.init_mesh()
TRS.init_systems()
TRS.setupSource(source_runtime=source_runtime, source_no_PML_cells=source_no_PML_cells)
TRS.set_dt(dt)
no_steps = int(N.ceil(runtime/dt))
TRS.runSteps(no_steps)
resses = {order:{dt_div:TRS.getResult()}}
drv_ts = N.array(TRS.inc_modeintg[::dt_div], N.float64)

#pickle.dump(resses, file('+brick_4_wg_tfsf_tmp.pickle', 'w'))
f_min, f_max = 1.05*f_c, 2*f_c
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

z_tf = z_measure - z_inc 
an_tf = N.array([an_calc.TF(f*2*N.pi, z_tf) for f in fr], N.complex128)

from NewCode.Consts import c0
line_types = ('k-', 'k-.', 'k--', 'k--o', 'k-.D', 'k-^', 'k--x')
"""
plot(fr*c0/1e9, 20*N.log10(an_tf), 'k', label='analytical |S21|')
plot(fr*c0/1e9, 20*N.log10(tf), line_types[1], label='3rd order |S21|')
plot(fr*c0/1e9, 20*N.log10(tf1), line_types[2], label='3rd order |S11|')
plot(fr*c0/1e9, 20*N.log10(tf), '--', label='2nd order |S21|')
plot(fr*c0/1e9, 20*N.log10(tf1), '--', label='2nd order |S11|')
"""
