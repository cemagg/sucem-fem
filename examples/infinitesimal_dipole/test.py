from __future__ import division

import numpy as N
import sys
import dolfin as dol
sys.path.append('../../')
sys.path.append('../../fenics_test/ABC/')
import test_data
import driver
reload(driver)
from FenicsCode.Consts import eps0, mu0, c0
from FenicsCode.Utilities.MeshIO import femmesh_2_dolfin_mesh
import dolfin as dol

order = 1
source_coord = N.array([0,0,0], N.float64) + 1e-5
source_pt = dol.Point(*source_coord)

parameters = dict(test_data.problem_data)
workspace = {}


lam = c0/test_data.problem_data['f']
# domain_size = N.array([4*lam]*3)
# max_edge_len = lam/6
# domain_subdivisions = N.array(N.ceil(domain_size/max_edge_len), N.uint)
# print domain_subdivisions
# #domain_subdivisions = N.array([21,21,21], N.uint)
# parameters.update(dict(domain_size=domain_size,
#                        domain_subdivisions=domain_subdivisions))
mesh_file = '../../workspace/sphere-r1m-6.femmesh'
mesh = femmesh_2_dolfin_mesh(mesh_file)
mesh.init()
workspace['mesh'] = mesh

parameters.update(dict(source_coord=source_coord,
                       order=order,
                       solver='direct'))
driver.run(parameters, workspace)
E_1 = driver.get_E_field(workspace, test_data.r_1)
#E_2 = driver.get_E_field(workspace, test_data.r_2)

V = workspace['V']
mesh = workspace['V'].mesh()
x = workspace['x']
A = workspace['A']
b = workspace['b']

u_abs = dol.Function(V)
u_re = dol.Function(V)
u_im = dol.Function(V)
u_abs.vector()[:] = N.abs(x)
u_re.vector()[:] = N.real(x)
u_im.vector()[:] = N.imag(x)
#dol.plot(u, interactive=True, axes=True)
#dol.plot(mesh, interactive=True, axes=True)

last = 18
ratio = N.abs(E_1[last] / test_data.E_1[last])

u_abs_recons = N.array([u_abs(r) for r in test_data.r_1[0:last+1]])
u_re_recons = N.array([u_re(r) for r in test_data.r_1[0:last+1]])
u_im_recons = N.array([u_im(r) for r in test_data.r_1[0:last+1]])



E_1_norm = E_1[0:last+1]/ratio


