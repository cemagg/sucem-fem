from __future__ import division

import numpy as N
import sys
sys.path.append('../../')
sys.path.append('../../fenics_test/ABC/')
import test_data
import driver
reload(driver)
from FenicsCode.Consts import eps0, mu0, c0

lam = c0/test_data.problem_data['f']
domain_size = N.array([2*lam]*3)
order = 1
max_edge_len = lam/15
source_coord = N.array([0,0,0], N.float64)
#domain_subdivisions = N.array(N.ceil(domain_size/max_edge_len), N.uint)
domain_subdivisions = N.array([5,5,5], N.uint)
parameters = dict(test_data.problem_data)
parameters.update(dict(domain_size=domain_size,
                       domain_subdivisions=domain_subdivisions,
                       source_coord=source_coord,
                       order=order))
workspace = {}
driver.run(parameters, workspace)
#E_1 = driver.get_E_field(workspace, r_1)
#E_2 = driver.get_E_field(workspace, r_2)

mesh = workspace['V'].mesh()
x = workspace['x']
x_lu = workspace['x_lu']
A = workspace['A']
ml = workspace['ml']
b = workspace['b']
print N.max(N.abs(1 - x/x_lu)), N.max(N.abs(x-x_lu))
