from __future__ import division

import numpy as N
import dolfin as dol
from dolfin import dot, curl, inner, dx, Point, PointSource

class CavityDims(object): pass
cdims = CavityDims()
cdims.a, cdims.b, cdims.c = 29,23,19

# Define mesh
mesh = dol.UnitCube(1,1,1)
mesh.coordinates()[:] *= [cdims.a,cdims.b,cdims.c]

# Define function space
order = 1
V = dol.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", order)

# Define basis and bilinear form
u = dol.TrialFunction(V)
v = dol.TestFunction(V)
source_coords = N.array([cdims.a/2, cdims.b/2, cdims.c/2], dtype=N.float64)
source_pt = Point(*source_coords)

m = inner(v, u)*dx                  # Mass form

M = dol.PETScMatrix()
dol.assemble(m, tensor=M, mesh=mesh)
b = dol.Vector(M.size(0))
b.array()[:] = 0

# ps0 = PointSource(V.sub(0), source_pt)
# ps1 = PointSource(V.sub(1), source_pt)
# ps2 = PointSource(V.sub(2), source_pt)

# ps0.apply(b)
# ps1.apply(b)
# ps2.apply(b)

dm = V.dofmap()
dofs = N.zeros(dm.max_cell_dimension(), dtype='I')
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
finite_element.evaluate_basis(el_basis_vals_all, source_coords, c)

xhat = N.zeros(4, N.float64)
