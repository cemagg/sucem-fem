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
from NewCode.GeomGen.Hybrid import ElementCloseDeadElements, WGSizeCalc, \
     MultiPinMinGeomlens
from NewCode.GeomGen.WGPinElements import CylindricalSection,\
     PointInCylindricalSection, PointInCylindricalSections
from RunInfos.WGPost.WR90HalfCylPosts import MultiCylPosts
from RunInfos.WGParms import WGParms
import NewCode.eMAGUSImport as eMAGUSImport
from NewCode.Analytical import WaveguidePhasor
from NewCode import Consts
try: max_dt_div_pow = int(sys.argv[1])
except IndexError: max_dt_div_pow = 8
try: FEKO_workspace = sys.argv[2]
except IndexError: FEKO_workspace = 'workspace'
try: output_dir = sys.argv[3]
except IndexError: output_dir = '.'
try: runfac = float(sys.argv[4])
except IndexError: runfac = 1.
    
eMAGUSImport.init(FEKO_workspace)
tet_mesh = eMAGUSImport.get_mesh()

print 'Tet_mesh elements: ', len(tet_mesh.elements)
tetmesh_basename, tetmesh_ext = os.path.splitext(tet_mesh.FemmeshFilename)
pin_name, h_div = re.compile(r'wg_(lech03fig.*?)_h(\d*).*').match(tetmesh_basename).groups()
h_div = int(h_div)
output_filename_template = output_dir + os.sep + 'hybmesh_%s_o%d_dt%e%s.pickle'
#pin_name = 'lech03fig7'

pin_parms = MultiCylPosts[pin_name]
WG = WGParms[pin_parms.GuideType]
a,b = WG.a, WG.b
an_calc = WaveguidePhasor.TE01(a)
h = a/h_div
f_c = an_calc.k_c/2/N.pi
lam_c = 1/f_c
drv_fun = Waveforms.get_gausspulse(fc=1.8*f_c, bw=1.2,tpr=-60)

r_ps = pin_parms.r_ps ; p_cs = pin_parms.p_cs ; phi_0s = pin_parms.phi_0s
r_es = N.zeros_like(r_ps) + 0.25*h
pinlen_info = MultiPinMinGeomlens(p_cs, r_ps, r_es, )
WGS = WGSizeCalc(h, pinlen_info.min_geomlen, pin_parms.min_buf_len)
z_measure = WGS.TF_len
pin_x_offs = a/2 if pin_parms.x_centered else 0.
pin_z_offs = z_measure/2 if pin_parms.z_centered else 0.
pin_offs = N.array([pin_x_offs, 0., pin_z_offs], N.float64)
pins = [CylindricalSection(p_c+pin_offs, r_p, phi_0)
        for p_c, r_p, phi_0 in izip(p_cs, r_ps, phi_0s)]
highres_in_pins = PointInCylindricalSections(pins, r_es)
lowres_in_pins = PointInCylindricalSections(
    pins, N.zeros_like(r_es)+N.sqrt(3*h**2)*0.51)

z_inc = 0
runtime = 62.5/f_c*runfac
source_runtime = runtime

no_PML_cells = 30
source_no_PML_cells = no_PML_cells

# dt_div = 1 corresponds to just-above stability dt for a/3 1st order element
base_dt = a/3*Consts.lumped_stability_factors[1]/1.1


orders = (2,)
#orders = (1,)
for order in orders:
    print 'order: ', order
    TRS = HybridWaveGuide2PortPEC(a,b,h,order,drv_fun, no_PML_cells)
    TRS.set_dims(z_inc, z_inc, z_measure)
    TRS.init_mesh()
    h_min = N.min(TRS.hex_mesh.gridStepSize)
    dt_div_min = int(N.ceil(base_dt/(h_min*Consts.lumped_stability_factors[order])))
    dt_divs = [dt_div_min] + [2**i for i in range(1,max_dt_div_pow+1)
                              if 2**i > dt_div_min]
    dt_divs = [32]
    hyb_geom = ElementCloseDeadElements(TRS.hex_mesh)
    hyb_geom.set_dead_node_selector(
        lowres_in_pins.in_pins, highres_in_pins.in_pins, [5,1,5])
    TRS.init_hyb_geom(hyb_geom)
    TRS.init_tetmesh(tet_mesh)
    TRS.init_hybrid_PEC(PEC_label=1001)
    TRS.init_systems()
    TRS.setupSource(source_runtime=source_runtime, source_no_PML_cells=source_no_PML_cells)
    
    for dt_div in dt_divs:
        dt = base_dt/dt_div 
        TRS.log_divisor = dt_div
        runfac_str = '_runfac%f' % runfac if runfac != 1 else ''
        output_filename = output_filename_template % (
            tetmesh_basename, order, dt, runfac_str)
        print dt_div, output_filename
        if os.path.exists(output_filename):
            print 'Output file "%s" exists, skipping' % output_filename
            continue
        TRS.set_dt(dt)
        no_steps = int(N.ceil(runtime/dt))
        TRS.runSteps(no_steps)
        res = TRS.getResult()
        res['drv_ts'] = N.array(TRS.inc_modeintg[::dt_div], N.float64)
        res['run_dt'] = dt
        res['base_dt'] = base_dt
        pickle.dump(res, file(output_filename, 'w'))
        


