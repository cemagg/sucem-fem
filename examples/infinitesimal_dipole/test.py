from __future__ import division

import numpy as N
import sys
import dolfin as dol
sys.path.append('../../')
import test_data
import driver
reload(driver)
from FenicsCode.Consts import eps0, mu0, c0
from FenicsCode.Utilities.MeshIO import femmesh_2_dolfin_mesh
import dolfin as dol

order = 2
source_coord = N.array([0,0,0.]) 
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
# mesh = dol.UnitCube(3,3,3)
# mesh.coordinates()[:] -= 0.5
mesh.coordinates()[:] *= lam
workspace['mesh'] = mesh

parameters.update(dict(source_coord=source_coord,
                       order=order,
                       solver='iterative'))
driver.run(parameters, workspace)
E_1 = driver.get_E_field(workspace, test_data.r_1)
last_E_1 = list(N.isnan(E_1[:,2])).index(True)
#E_2 = driver.get_E_field(workspace, test_data.r_2)

V = workspace['V']
mesh = workspace['V'].mesh()
x = workspace['x']
A = workspace['A']
b = workspace['b']

M = workspace['M']
S = workspace['S']
S_0 = workspace['S_0']


from pylab import *
r1 = test_data.r_1[1:last_E_1+1]/lam
x1 = r1[:,0]
E_1_ana = test_data.E_1[1:last_E_1+1]
E_1_num = E_1[1:last_E_1+1]
figure(1)
plot(x1, N.abs(E_1_num[:,0]), '-g', label='x_num')
plot(x1, N.abs(E_1_num[:,1]), '-b', label='y_num')
plot(x1, N.abs(E_1_num[:,2]), '-r', label='z_num')
plot(x1, N.abs(E_1_ana[:,2]), '--r', label='z_ana')
ylabel('E-field Magnitude')
xlabel('Distance (wavelengths)')
legend(loc='best')
show()
