from __future__ import division
import numpy as N

from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.GeomGen.Hybrid import SimpleSurroundedRect
from NewCode.Analytical import WaveguidePhasor



geom_filename = '+WG_tst.geo'

a,b = 22.86e-3, 10.16e-3
an_calc = WaveguidePhasor.TE01(a)

f_c = an_calc.k_c/2/N.pi
lam_c = 1/f_c
h = (lam_c/2)/5.
z_inc = 0
z_measure = (lam_c/2)*7/5
c = z_measure 
hex_mesh = BrickMesh.Mesh(BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c,h))

pin_z_cent = a*(1/2 + 1/5)
pin_x_cent = a/5
pin_d = a/5

wiggle_room = N.array([h/2, 0, h/2], N.float64)

bot_left = N.array([a/2-pin_d/2+pin_x_cent, 0, pin_z_cent-pin_d/2],
                   N.float64) - wiggle_room
top_right = N.array([a/2+pin_d/2+pin_x_cent, b, pin_z_cent+pin_d/2],
                    N.float64) + wiggle_room


GG = SimpleSurroundedRect(hex_mesh)
GG.set_tet_range(bot_left, top_right)
GG.init_dead_background_element_set()
GG.init_geom(include_meshbdry_els=False)
GG.write_geom(file(geom_filename, 'w'), write_vols=False)

