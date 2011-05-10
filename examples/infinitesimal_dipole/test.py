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

lam = c0/test_data.problem_data['f']
domain_size = N.array([4*lam]*3)
order = 2
max_edge_len = lam/3
source_coord = N.array([0,0,0], N.float64)
domain_subdivisions = N.array(N.ceil(domain_size/max_edge_len), N.uint)
print domain_subdivisions
#domain_subdivisions = N.array([21,21,21], N.uint)
parameters = dict(test_data.problem_data)
parameters.update(dict(domain_size=domain_size,
                       domain_subdivisions=domain_subdivisions,
                       source_coord=source_coord,
                       order=order))
workspace = {}
driver.run(parameters, workspace)
E_1 = driver.get_E_field(workspace, test_data.r_1)
#E_2 = driver.get_E_field(workspace, test_data.r_2)

V = workspace['V']
mesh = workspace['V'].mesh()
x = workspace['x']
A = workspace['A']
b = workspace['b']

u = dol.Function(V)
u.vector()[:] = x
#dol.plot(u, interactive=True)
#dol.plot(mesh, interactive=True, axes=True)

last = 18
ratio = N.abs(E_1[last] / test_data.E_1[last])

E_1_norm = E_1[0:last+1]/ratio
