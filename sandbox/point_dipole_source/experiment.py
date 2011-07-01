# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import numpy as N
import dolfin as dol
from dolfin import dot, curl, inner, dx, Point, PointSource

class CavityDims(object): pass
cdims = CavityDims()
cdims.a, cdims.b, cdims.c = 1,1,1

# Define mesh
mesh = dol.UnitCube(1,1,1)
mesh.coordinates()[:] *= [cdims.a,cdims.b,cdims.c]

# Define function space
order = 1
V = dol.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", order)

# Define basis and bilinear form
u = dol.TrialFunction(V)
v = dol.TestFunction(V)
source_coords = N.array([cdims.a/2, cdims.b/2, cdims.c/2], N.float64)
source_pt = Point(*source_coords)
source_value = N.array([0,0,1], N.float64)

m = inner(v, u)*dx                  # Mass form

M = dol.PETScMatrix()
dol.assemble(m, tensor=M, mesh=mesh)
b = dol.Vector(M.size(0))
b.array()[:] = 10

dm = V.dofmap()
dofs = N.zeros(dm.max_cell_dimension(), dtype=N.uint32)
cell_index = mesh.any_intersected_entity(source_pt)
#cell_index = 
c = dol.Cell(mesh, cell_index)
# Check that the source point is in this element    
assert(c.intersects_exactly(source_pt)) 

dm.tabulate_dofs(dofs,  c)
finite_element = V.dolfin_element()
no_basis_fns = finite_element.space_dimension()
# Vector valued elements have rank of 1
assert(finite_element.value_rank() == 1)
# Vector valued elements have only one rank (i.e. 0) along which
# dimensions are defined. This is the dimension that the basis
# function value vector is. Since we have 3D Nedelec elements here
# this should be 3
bf_value_dimension = finite_element.value_dimension(0)
el_basis_vals = N.zeros((no_basis_fns, bf_value_dimension), dtype=N.float64)
finite_element.evaluate_basis_all(el_basis_vals, source_coords, c)
rhs_contribs = N.sum(el_basis_vals*source_value, axis=1)
b[dofs] = rhs_contribs


## Below don't work due to incomplete python interface
# xhat = N.zeros(4, N.float64)
# finite_element.map_to_reference_cell(xhat, source_coords, c)
##
