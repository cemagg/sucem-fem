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
