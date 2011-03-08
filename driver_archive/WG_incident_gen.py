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
from scipy.io import numpyio
#
# Local Imports
#
import NewCode
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import Waveforms
from NewCode.Feeds.WaveGuideTFSFBootstrap import WaveGuideTFSFBootstrap
from NewCode.IOUtils import read_filelog

from analytical_WG_driver import WGT

dt_divisors = {1:(1,128),               
               2:(3,128),
               3:(4,128),
               4:(128,9)}
omega_c = N.pi                          # Lowest freq cutoff mode frequency for
                                        # WG with largest cross-section of 1m    
a,b = 1., 0.25, 
output_filebasename = '+_WG_inc'
h = 1/8.
order = 3
dt = 1/16/dt_divisors[order][0]
runtime = 10
no_PML_cells = 15
drv_fun = WGT.discrete_drv_fn

output_filename = output_filebasename + "%sX%s_h%s_o_%s_dt%s_pml%s" %(
    str(a), str(b), str(1/h), str(order), str(1/dt), str(no_PML_cells))
E_logfilename = output_filename+'_E.numpyraw'
B_logfilename = output_filename+'_B.numpyraw'
E_logfile = file(E_logfilename, 'wb')
B_logfile = file(B_logfilename, 'wb')


print 'order: ', order
WGB = WaveGuideTFSFBootstrap(a,b,h,order, drv_fun, no_PML_cells=no_PML_cells)
WGB.init_mesh()
WGB.init_systems()
WGB.setupSource()
WGB.set_dt(dt)
WGB.init_steppers(runtime, yield_zero=False)

type_char = 'd'
E_doflen = len(WGB.log_dofnos.E)
B_doflen = len(WGB.log_dofnos.B)
E_logfile.write(type_char + ' ' + str(E_doflen) + ' \n')
B_logfile.write(type_char + ' ' + str(B_doflen) + ' \n')

mode_intg = []
while True:
    try: res = WGB.next_drive()
    except StopIteration: break
    numpyio.fwrite(E_logfile, E_doflen, res.E_incdofs)
    numpyio.fwrite(B_logfile, B_doflen, res.B_incdofs)
    mode_intg.append(res.mode_weighted)

pickle.dump(Struct(mode_intg=mode_intg), file(output_filename + '.pickle', 'w'))
E_logfile.close()
B_logfile.close()
            
E_logged = read_filelog(file(E_logfilename, 'rb'))
B_logged = read_filelog(file(B_logfilename, 'rb'))
