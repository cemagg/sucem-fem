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
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.Runners.HybridWaveGuide2Port import HybridWaveGuide2PortPEC
from NewCode.GeomGen.Hybrid import ElementCloseDeadElements, WGSizeCalc, \
     MultiPinMinGeomlens
from NewCode.GeomGen.WGPinElements import CylindricalSection,\
     PointInCylindricalSection, GenPinsGeo, PointInCylindricalSections
from RunInfos.WGParms import WGParms
from RunInfos.WGPost.WR90HalfCylPosts import MultiCylPosts
import NewCode.eMAGUSImport as eMAGUSImport
from NewCode.Analytical import WaveguidePhasor
from NewCode import Consts

h_div = 8.
geom_filename = 'mesh_geoms/wg_lech03fig4r4.41_h%d.geo' % int(h_div)
#geom_filename = 'mesh_geoms/wg_lech03fig7_h%d.geo' % int(h_div)
WG = WGParms['WR90']
#pin_name = 'lech03fig7'
pin_name = 'lech03fig4r4.41'
pin_parms = MultiCylPosts[pin_name]
a,b = WG.a, WG.b
h = a/h_div

r_ps = pin_parms.r_ps ; p_cs = pin_parms.p_cs ; phi_0s = pin_parms.phi_0s
r_es = N.zeros_like(r_ps) + 0.25*h
pinlen_info = MultiPinMinGeomlens(p_cs, r_ps, r_es)
WGS = WGSizeCalc(h, pinlen_info.min_geomlen, pin_parms.min_buf_len)
z_measure = WGS.TF_len
c = z_measure
hex_mesh = BrickMesh.Mesh(BrickMeshGen.make_rect_cavity_brick_listmesh(
    a,b,c,h))
print 'Hex Mesh elements: ', len(hex_mesh.elements)

pin_x_offs = a/2 if pin_parms.x_centered else 0.
pin_z_offs = z_measure/2 if pin_parms.z_centered else 0.
pin_offs = N.array([pin_x_offs, 0., pin_z_offs], N.float64)
pins = [CylindricalSection(p_c+pin_offs, r_p, phi_0)
        for p_c, r_p, phi_0 in izip(p_cs, r_ps, phi_0s)]
highres_in_pins = PointInCylindricalSections(pins, r_es)
lowres_in_pins = PointInCylindricalSections(
    pins, N.zeros_like(r_es)+N.sqrt(3*h**2)*0.51)
pins_geo = GenPinsGeo(pins, b)
pins_geo.calc_nodes()
pins_geo.calc_edges()
pins_geo.calc_segs()


GG = ElementCloseDeadElements(hex_mesh)
GG.set_dead_node_selector(lowres_in_pins.in_pins, highres_in_pins.in_pins, [5,1,5])
GG.init_dead_background_element_set()
GG.init_geom(include_meshbdry_els=False)
GG.add_extra_geom_nodes(pins_geo.nodes)
GG.add_extra_geom_edges(pins_geo.edges)
GG.add_extra_geom_circsegs(pins_geo.segs)
GG.write_geom(file(geom_filename, 'w'), h_geo=.75*h, write_vols=False)
