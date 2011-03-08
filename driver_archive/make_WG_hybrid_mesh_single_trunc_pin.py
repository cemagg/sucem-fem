from __future__ import division
import numpy as N

from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.GeomGen.WGPinElements import CylindricalSection,\
     PointInCylindricalSection, GenPinGeo
from NewCode.GeomGen.Hybrid import NodeInDeadElements, WGSizeCalc
from RunInfos.WGPost.WR90HalfCylPosts import SingleCylPosts
from RunInfos.WGParms import WGParms


geom_filename = '+WG_tst.geo'

pin_parms = SingleCylPosts['hcyl_pin_1']
WG = WGParms[pin_parms.GuideType]
a,b = WG.a, WG.b

h_div = 64.
h = a/h_div
z_inc = 0
r_e = 0.25*h
pin_min_geomlen = 2*r_e + 2*pin_parms.r
min_buflen = pin_parms.min_buf_len
WGS = WGSizeCalc(h, pin_min_geomlen, min_buflen)
z_measure = WGS.TF_len
c = z_measure 
hex_mesh = BrickMesh.Mesh(BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c,h))
print 'Hex Mesh elements: ', len(hex_mesh.elements)
pin_cent = [a/2-pin_parms.d, 0, z_measure/2]
pin = CylindricalSection(pin_cent, pin_parms.r, phi_0=pin_parms.phi)
in_pin = PointInCylindricalSection(pin, r_e)

pin_geo = GenPinGeo(pin, b)
pin_geo.calc_nodes()
pin_geo.calc_edges()
pin_geo.calc_segs()

GG = NodeInDeadElements(hex_mesh)
GG.set_dead_node_selector(in_pin.in_pin)
GG.init_dead_background_element_set()
GG.init_geom(include_meshbdry_els=False)
GG.add_extra_geom_nodes(pin_geo.nodes)
GG.add_extra_geom_edges(pin_geo.edges)
GG.add_extra_geom_circsegs(pin_geo.segs)
GG.write_geom(file(geom_filename, 'w'), h_geo=.75*h, write_vols=False)

