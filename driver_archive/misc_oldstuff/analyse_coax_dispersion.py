"""
Generates k-w diagram assuming a reference phase of 0 at z=0 and a measured
field point at z = p0 m
"""
from __future__ import division
from itertools import izip
import os, sys
import pickle
import random
import numpy as N
import scipy
from NewCode import Waveforms
from NewCode.Utilities import Struct

p0 = 5.                                 # z-coord of measurement point

drv_fun = Waveforms.get_d_gaussian(fc=0.5, tpr=-50)
dt = 0.01
import analyse_fftstuff
class AnalyseThingyBoxcar(analyse_fftstuff.AnalyseThingy):
    windowFun = staticmethod(scipy.signal.boxcar)
#rsvals should b
file_basename = '+%s_%s_rsvals.pickle'
tests = [
    ('newmark_coax', '0.3', 'Wave Eqn', 'k-'),
    ('coupleda_coax', '0.3', "EBHD Maxwell's", 'k-.'),
    #('coupleda_coax', '0.15'),
    #('coupleda_coax', '0.1'),
    ('coupledb_coax', '0.3', "EB Maxwell's", 'k--'),
    #('coupledb_coax', '0.15'),
    #('coupledb_coax', '0.1'),
]

xmin, xmax = 50, 250
test_results = dict()
for testname, meshdensity, label, style in tests:
    rsvals = pickle.load(file(file_basename % (testname, meshdensity)))
    at = AnalyseThingyBoxcar([rsvals[:,200,0]], dt, drv_fun)
    at.xmin = xmin
    at.xmax = xmax
    TF = at.dfts_series[0]/at.fs
    k_calc =  -N.unwrap(N.angle(TF[0:at.xmax]))[at.xmin:at.xmax]/p0
    test_results[(testname, meshdensity)] = Struct(
        rsvals=rsvals, at=at, TF=TF, k_calc=k_calc)

freq_range = N.arange(at.xmin,at.xmax)*at.df

# k is in radians/meter
k_analytical = freq_range*2*N.pi

for testname, meshdensity, label, style in tests:
    plot(freq_range*N.pi*2, test_results[(testname, meshdensity)].k_calc,
         style, label=label)

plot(freq_range*N.pi*2, k_analytical, 'k:', label='analytical')
xlabel('Frequency (rad/s)')
ylabel('Wavenumber (rad/m)')
legend(loc=2)
show()
