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
from NewCode.Utilities import Struct, partial, close_to_point, in_box_vec, in_box
from NewCode import Waveforms
from NewCode.Runners.WaveGuide2Port import WaveGuide2PortPECelsBuff
from NewCode.Analytical import WaveguidePhasor
from RunInfos.WGPost.WR90HalfCylPosts import SingleCylPosts
from RunInfos.WGParms import WGParms
from NewCode.GeomGen.Hybrid import DeadElements
from NewCode.GeomGen.WGPinElements import CylindricalSection,\
     PointInCylindricalSection, GenPinGeo
from NewCode import Consts

pin_name = 'hcyl_pin_1'
output_dir = 'hcyl_pin_output'
# brick_xy... x = b -> minimum buffer space, y = n -> all nodes inside 
# brick_b   ... buffer space, centroid inside
# brick   ... no buffer space, centroid inside
output_filename_template = output_dir+os.sep+'brick_bn_wg_hcyl_pin_1_h%d_o%d_dt_%f.pickle'
#output_filename_template = output_dir+os.sep+'brick_b_wg_hcyl_pin_1_h%d_o%d_dt_%f.pickle'
nodesel_type = re.match(r'.*brick_.*?(n?)_wg.*', output_filename_template).groups()[0]
#h_divs = [8, 16, 32, 48, 64]
#h_divs = [8, 11, 16, 21, 32, 48, 64, 95, 128]
h_divs = [128]
order = 1
max_dt_div_pow = 1
pin_parms = SingleCylPosts[pin_name]
WG = WGParms[pin_parms.GuideType]
a,b = WG.a, WG.b
an_calc = WaveguidePhasor.TE01(a)

f_c = an_calc.k_c/2/N.pi
lam_c = 1/f_c
drv_fun = Waveforms.get_gausspulse(fc=1.8*f_c, bw=1.2,tpr=-60)
# dt_div = 1 corresponds to just-above stability dt for a/3 1st order element
base_dt = a/3*Consts.lumped_stability_factors[1]/1.1

z_inc = 0
drv_fun = Waveforms.get_gausspulse(fc=1.5*f_c, bw=1.1,tpr=-60)

runtime = 62.5/f_c
source_runtime = runtime
no_PML_cells = 30
source_no_PML_cells = no_PML_cells
min_bufflen = a/16*30

for h_div in h_divs:
    h = a/h_div
    z_measure = int(N.ceil(((lam_c/2)*5/4)/h))*h
    pin_cent = [a/2-pin_parms.d, 0, z_measure/2]
    pin = CylindricalSection(pin_cent, pin_parms.r, phi_0=pin_parms.phi)
    in_pin = PointInCylindricalSection(pin)
    z_PML_goal = z_measure
    TRS = WaveGuide2PortPECelsBuff(a,b,h,order,drv_fun, no_PML_cells)
    buffspace = 0.
    if h*no_PML_cells < min_bufflen:
        buffspace += min_bufflen - h*no_PML_cells
    print 'buffspace: %f' % buffspace
    TRS.set_dims(z_inc, z_inc, z_measure, buffspace=buffspace)
    TRS.init_mesh()
    if nodesel_type == '':
        print "Centroid in pin selected"
        pin_els = set(i for i, inel in enumerate(in_pin.in_pin(elc) for elc in N.average(
            TRS.hex_mesh.elements[:].nodeCoords, axis=1)) if inel)
    elif nodesel_type == 'n':
        print "Any node in pin selected"
        pin_els = set(i for i, inel in enumerate(
            N.all([in_pin.in_pin(co) for co in elc])
            for elc in TRS.hex_mesh.elements[:].nodeCoords)
                      if inel)
    else: raise Exception("Unknown pin-element selection method")
        
#     pin_geo = GenPinGeo(pin, b)
#     pin_geo.calc_nodes()
#     pin_geo.calc_edges()
#     pin_geo.calc_segs()

#     GG = DeadElements(TRS.hex_mesh)
#     GG.dead_background_elements = pin_els
#     GG.init_geom(include_meshbdry_els=True)
#     GG.add_extra_geom_nodes(pin_geo.nodes)
#     GG.add_extra_geom_edges(pin_geo.edges)
#     GG.add_extra_geom_circsegs(pin_geo.segs)
#     GG.write_geom(file('+WG_brickpin_n_h%d.geo' %h_div, 'w'), h_geo=.75*h, write_vols=False)
    h_min = N.min(TRS.hex_mesh.gridStepSize)
    dt_div_min = int(N.ceil(base_dt/(h_min*Consts.lumped_stability_factors[order])))
    dt_divs = [dt_div_min] + [2**i for i in range(1,max_dt_div_pow+1)
                              if 2**i > dt_div_min]
    TRS.set_PECels(pin_els)
    TRS.init_systems()
    TRS.setupSource(source_runtime=source_runtime, source_no_PML_cells=source_no_PML_cells)
    for dt_div in dt_divs:
        dt = base_dt/dt_div
        output_filename = output_filename_template % (h_div, order, dt)
        print dt_div, output_filename
#         if os.path.exists(output_filename):
#             print 'Output file "%s" exists, skipping' % output_filename
#             continue
        TRS.log_divisor = dt_div
        TRS.set_dt(dt)
        no_steps = int(N.ceil(runtime/dt))
        TRS.runSteps(no_steps)
        res = TRS.getResult()
        res['drv_ts'] = N.array(TRS.inc_modeintg[::dt_div], N.float64)
        res['run_dt'] = dt
        res['base_dt'] = base_dt
        #pickle.dump(res, file(output_filename, 'w'))

