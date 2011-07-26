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

import sys
sys.path.append('../../')
import numpy as N
import os
import dolfin as dol
from dolfin import dot, curl, inner, dx
from scipy.sparse.linalg.eigen.arpack import speigs

from sucemfem.Consts import eps0, mu0, c0

import postproc_eigres
from material_properties import MaterialProperties

# Test for PETSc and SLEPc
if not dol.has_la_backend("PETSc"):
    print "DOLFIN has not been configured with PETSc. Exiting."
    exit()

if not dol.has_slepc():
    print "DOLFIN has not been configured with SLEPc. Exiting."
    exit()

# class CavityDims(object): pass
# cdims = CavityDims()
# cdims.a, cdims.b, cdims.c = 29,23,19


# Define mesh
# mesh = dol.UnitCube(1,1,1)
# mesh.coordinates()[:] *= [cdims.a,cdims.b,cdims.c]
#mesh_file = 'lee_mittra92_fig6b.xml'
#mesh_file = 'lee_mittra92_fig6c.xml'
mesh_file = '../examples/albani_bernardi74/mesh/albani_bernardi74_fig2VII.xml'
materials_mesh_file = "%s_physical_region%s" % (os.path.splitext(mesh_file))
mesh = dol.Mesh(mesh_file)
materials = {1000:MaterialProperties(eps_r=16/eps0),
             1001:MaterialProperties(eps_r=1/eps0)}
eps_vals = dict((k, mat.get_eps()) for k,mat in materials.items())

material_mesh_func = dol.MeshFunction('uint', mesh, materials_mesh_file)
# Discontinuous function to represent material_properties
material_func_space = dol.FunctionSpace(mesh, 'DG', 0)
eps = dol.Function(material_func_space)
eps.vector()[:] = N.array([eps_vals[int(i)] for i in material_mesh_func.array()])


# Define function space
order = 3
V = dol.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", order)

# Define basis and bilinear form
u = dol.TrialFunction(V)
v = dol.TestFunction(V)

m = eps*inner(v, u)*dx                  # Mass form
# m = inner(v, u)*dx                  # Mass form
s = dot(curl(v), curl(u))*dx            # Stiffness form
#s = (1/mu0)*dot(curl(v), curl(u))*dx            # Stiffness form

def boundary(x, on_boundary):
    return on_boundary
zero = dol.Expression(("0.0", "0.0", "0.0"), degree=1)
bc = dol.DirichletBC(V, zero, boundary)



# Assemble smass form
M = dol.PETScMatrix()
S = dol.PETScMatrix()
dol.assemble(m, tensor=M, mesh=mesh)
bc.apply(M)
dol.assemble(s, tensor=S, mesh=mesh)
bc.apply(S)
#M, S = dol.assemble_system(m, s, bc)


sigma = 1.5
smat = S - sigma*M
#lu = dol.LUSolver(S - sigma*M)
lu = dol.LUSolver(smat)
lu.parameters["reuse_factorization"] = True
lu.parameters["report"] = False
bb = dol.Vector(M.size(0))
xx = dol.Vector(M.size(0))
def sigma_solve(b):
    bb[:] = b
    lu.solve(xx, bb)
    return xx[:]
M_matvec = lambda x: M*x

arpack_eigs,arpack_v = speigs.ARPACK_gen_eigs(M_matvec, sigma_solve, M.size(0), sigma, 51, ncv=91)

# Create eigensolver
# esolver = dol.SLEPcEigenSolver(S,M)
# esolver.parameters["spectrum"] = "smallest real"
# esolver.parameters["solver"] = "arnoldi"
# esolver.parameters["spectral_shift"] = sigma
# esolver.parameters["spectral_transform"] = "shift-and-invert"
# esolver.solve()

# eigs = [esolver.get_eigenvalue(i)[0] for i in range(4)]
# filtered_eigs = []
# for i in range(S.size(0)):
#     ev_r, ev_i = esolver.get_eigenvalue(i)
#     if ev_r > 0.001:
#         filtered_eigs.append(ev_r)

res = N.array(sorted(arpack_eigs)[0:10])
print N.sqrt(res)
print c0*N.sqrt(res)/2/N.pi/1e6
#errs = postproc_eigres.calc_errs(res)
#postproc_eigres.print_errs(errs)


# class TestSub(dol.SubDomain):
#     retvals = [True, True, True, True, True]
#     i = 0
#     def inside(self, x, on_boundary):
#         #rv = self.retvals[self.i%len(self.retvals)]
#         rv = x[0] > -0.49
#         self.i += 1
#         return rv
# ts = TestSub()
# mf = dol.MeshFunction('uint', mesh, mesh.topology().dim())
# ts.mark(mf,1111)
# print mf.array()
