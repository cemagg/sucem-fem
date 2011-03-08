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
import re
#
# Local Imports
#
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import Waveforms
from NewCode.Runners.HybridWaveGuide2Port import HybridWaveGuide2PortPEC
from NewCode.GeomGen.Hybrid import NodeInDeadElements, WGSizeCalc
from NewCode.GeomGen.WGPinElements import CylindricalSection,\
     PointInCylindricalSection
from RunInfos.WGPost.WR90HalfCylPosts import SingleCylPosts
from RunInfos.WGParms import WGParms
import NewCode.eMAGUSImport as eMAGUSImport
from NewCode.Analytical import WaveguidePhasor
from NewCode import Consts

try:
    max_dt_div_pow = int(sys.argv[1])
except IndexError:
    max_dt_div_pow = 8
try:
    FEKO_workspace = sys.argv[2]
except IndexError:
    FEKO_workspace = 'workspace'
try:
    output_dir = sys.argv[3]
except IndexError:
    output_dir = '.'
    
eMAGUSImport.init(FEKO_workspace)
tet_mesh = eMAGUSImport.get_mesh()

print 'Tet_mesh elements: ', len(tet_mesh.elements)
tetmesh_basename, tetmesh_ext = os.path.splitext(tet_mesh.FemmeshFilename)
pinno, h_div = re.compile(r'wg_hcyl_pin_(\d*).*_h(\d*).*').match(tetmesh_basename).groups()
h_div = int(h_div)
pin_name = 'hcyl_pin_' + pinno
output_filename_template = output_dir + os.sep + 'time_hybmesh_%s_o%d_dt%f.pickle'

pin_parms = SingleCylPosts[pin_name]
WG = WGParms[pin_parms.GuideType]
a,b = WG.a, WG.b
an_calc = WaveguidePhasor.TE01(a)
h = a/h_div
f_c = an_calc.k_c/2/N.pi
lam_c = 1/f_c
drv_fun = Waveforms.get_gausspulse(fc=1.8*f_c, bw=1.2,tpr=-60)

r_e = 0.25*h
pin_min_geomlen = 2*r_e + 2*pin_parms.r
min_buflen = pin_parms.min_buf_len
WGS = WGSizeCalc(h, pin_min_geomlen, min_buflen)
z_measure = WGS.TF_len

pin_cent = [a/2-pin_parms.d, 0, z_measure/2]
pin = CylindricalSection(pin_cent, pin_parms.r, phi_0=pin_parms.phi)
in_pin = PointInCylindricalSection(pin, r_e)

l_sc = h
z_inc = 0

runtime = 62.5/f_c
source_runtime = runtime

no_PML_cells = 30
source_no_PML_cells = no_PML_cells
min_bufflen = a/16*30

# dt_div = 1 corresponds to just-above stability dt for a/3 1st order element
base_dt = a/3*Consts.lumped_stability_factors[1]/1.1


#orders = (1,3,2)
#orders = (3,)
orders = (3,)
for order in orders:
    print 'order: ', order
    TRS = HybridWaveGuide2PortPEC(a,b,h,order,drv_fun, no_PML_cells)
    buffspace = 0.
    if h*no_PML_cells < min_bufflen:
        buffspace += min_bufflen - h*no_PML_cells
    print 'buffspace: %f' % buffspace
    #TRS.set_dims(z_inc, z_inc, z_measure)
    TRS.set_dims(z_inc, z_inc, z_measure, buffspace=buffspace)
    TRS.init_mesh()
    h_min = N.min(TRS.hex_mesh.gridStepSize)
    dt_div_min = int(N.ceil(base_dt/(h_min*Consts.lumped_stability_factors[order])))
    #     dt_divs = [dt_div_min] + [2**i for i in range(1,max_dt_div_pow+1)
    #                               if 2**i > dt_div_min]
    #dt_divs = [dt_div_min]
    dt_divs = [32]
    hyb_geom = NodeInDeadElements(TRS.hex_mesh)
    hyb_geom.set_dead_node_selector(in_pin.in_pin)
    TRS.init_hyb_geom(hyb_geom)
    TRS.init_tetmesh(tet_mesh)
    TRS.init_hybrid_PEC(PEC_label=1001)
    TRS.init_systems()
    TRS.setupSource(source_runtime=source_runtime, source_no_PML_cells=source_no_PML_cells)
    for dt_div in dt_divs:
        dt = base_dt/dt_div 
        TRS.log_divisor = dt_div
        output_filename = output_filename_template % (tetmesh_basename, order, dt)
        TRS.hybrid_system.writestuff_n = 1200
        TRS.hybrid_system.writefile = output_filename[0:-7]
        TRS.hybrid_system.writestuff_exit = True
        print dt_div, output_filename
#         if os.path.exists(output_filename):
#             print 'Output file "%s" exists, skipping' % output_filename
#             continue
        TRS.set_dt(dt)
        if max_dt_div_pow == 0: runtime = dt
        no_steps = int(N.ceil(runtime/dt))
        TRS.runSteps(no_steps)
        res = TRS.getResult()
        res['drv_ts'] = N.array(TRS.inc_modeintg[::dt_div], N.float64)
        res['run_dt'] = dt
        res['base_dt'] = base_dt
        if max_dt_div_pow != 0: pickle.dump(res, file(output_filename, 'w'))
        


# # z_tf = z_measure - z_inc 
# # an_tf = N.array([an_calc.TF(f*2*N.pi, z_tf) for f in fr], N.complex128)
# import FEKO_WR90_hcyl_pin123 as FD
"""
figure(1); plot(FD.fr, FD.pin1_S11, label='FEKO')
figure(1); plot(fr, 20*N.log10(tf1), label='o%d, dt%f'% (order, dt))
"""
