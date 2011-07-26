## Copyright (C) 2011 Stellenbosch University
##
## This file is part of SUCEM.
##
## SUCEM is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SUCEM is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with SUCEM. If not, see <http://www.gnu.org/licenses/>. 
##
## Contact: cemagga@gmail.com 
from __future__ import division

"""Generate data used as reference values for the fillament source test"""

import dolfin
import numpy as np
import pickle
import sys
sys.path.insert(0, '../../../')

from sucemfem.Utilities.MeshGenerators import get_centred_cube
import sucemfem.Utilities.Optimization
from sucemfem.Consts import c0
from sucemfem.Sources.fillament_current_source import FillamentCurrentSource



sucemfem.Utilities.Optimization.set_dolfin_optimisation(True)


## Postprocessing requests
theta_deg = np.linspace(0, 180, 181)
no_ff_pts = len(theta_deg)
phi_deg = np.zeros(no_ff_pts)


## Problem parameters
freq =  1.0e+9                          # Frequency
lam = c0/freq
l = lam/4                               # Dipole length
# l = 2.3*lam                               # Dipole length
I = 1.0                                 # Dipole current
source_direction = np.array([0,0,1.])    # Source orientation
source_centre = np.array([0,0,0.])        # Position of the source
source_endpoints =  np.array(
    [-source_direction*l/2, source_direction*l/2]) + source_centre

## Discretisation settings
order = 2
domain_size = np.array([lam]*3)/2
max_edge_len = lam/3
mesh = get_centred_cube(domain_size, max_edge_len)
## Implementation
V = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", order)
dipole_source = FillamentCurrentSource()
dipole_source.no_integration_points = 20
dipole_source.set_function_space(V)
dipole_source.set_source_endpoints(source_endpoints)
dipole_source.set_value(I)
dofnos, rhs_contribs = dipole_source.get_contribution()
sortind = np.argsort(dofnos)
dofnos = dofnos[sortind].copy()
rhs_contribs = rhs_contribs[sortind].copy()

pickle.dump(dict(order=order,I=I, source_endpoints=source_endpoints,
                 domain_size=domain_size, max_edge_len=max_edge_len,
                 no_integration_points=dipole_source.no_integration_points,
                 dofnos=dofnos, rhs_contribs=rhs_contribs),
            open('fillament_source_test_data.pickle', 'w'))
