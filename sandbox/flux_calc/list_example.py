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
# Authors:
# Neilen Marais <nmarais@gmail.com>

from __future__ import division

import numpy as np
import dolfin

def as_dolfin_vector(a):
    """Convert array to a dolfin Vector() instance"""
    assert len(a.shape) == 1            # 1D vectors please
    v = dolfin.Vector(len(a))
    v.set_local(np.require(a, requirements=['C',]))
    return v


mesh = dolfin.UnitCube(1,1,1)
V = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", 2)
no = V.dofmap().global_dimension()
x = np.ones(no) 
u = dolfin.Function(V, as_dolfin_vector(x))
vec = as_dolfin_vector(x)
uu = dolfin.Function(V, vec)
