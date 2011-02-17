from __future__ import division
import numpy as N

from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.GeomGen.Hybrid import SimpleSurroundedRect

a,b,c = 1., 0.25, 1.
h = 1/5.
hex_mesh = BrickMesh.Mesh(BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c,h))

geom_filename = '+WG_tst.geo'

x_c = 1/3
z_c = 5/8
d=1/5
bot_left = N.array([x_c-d/2, 0, z_c-d/2], N.float64)
top_right = N.array([x_c+d/2, b, z_c+d/2], N.float64)

GG = SimpleSurroundedRect(hex_mesh)
GG.set_tet_range(bot_left, top_right)
GG.init_dead_background_element_set()
GG.init_geom()
GG.write_geom(file(geom_filename, 'w'))

