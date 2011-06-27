from __future__ import division

import sys
sys.path.append('../../')
import numpy as N
import os
import dolfin as dol
from dolfin import dot, cross, curl, inner, dx, ds
import scipy.sparse

from FenicsCode.Consts import eps0, mu0, c0
import point_source
reload(point_source)

# Define material and frequency
eps_r = 1;
mu_r = 1;
freq = 1e9;
source_coord = [0.5,0.5,0.5]
source_value = N.array([0,0,1.])

# Define mesh
mesh = dol.UnitCube(11,11,11)

# Define function space
order = 1
V = dol.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", order)

# Define basis and bilinear form
k_0 = 2*N.pi*freq/c0

u = dol.TrialFunction(V)
v = dol.TestFunction(V)

m = eps_r*inner(v, u)*dx                  # Mass form
s = (1/mu_r)*dot(curl(v), curl(u))*dx            # Stiffness form

n = V.cell().n
s_0 = inner ( cross ( n, v), cross ( n, u ) )*ds

def boundary(x, on_boundary):
    return on_boundary

# Assemble forms
M = dol.uBLASSparseMatrix()
S = dol.uBLASSparseMatrix()
S_0 = dol.uBLASSparseMatrix()

dol.assemble(m, tensor=M, mesh=mesh)
dol.assemble(s, tensor=S, mesh=mesh)
dol.assemble(s_0, tensor=S_0, mesh=mesh)

from scipy.sparse import csr_matrix
from pyamg import smoothed_aggregation_solver
from numpy import intc

def to_scipy_csr(mat, dtype=None, imagify=False):
    (row,col,data) = mat.data()   # get sparse data
    col = intc(col)
    row = intc(row)
    n = mat.size(0)
    if imagify: data = data*1j
    Asp = csr_matrix( (data,col,row), shape=(n,n), dtype=dtype)
    return Asp

Msp = to_scipy_csr(M)
Ssp = to_scipy_csr(S)
S_0sp = to_scipy_csr(S_0, dtype=N.complex128, imagify=True)

A = Ssp - k_0**2*Msp + k_0*S_0sp 
b = N.zeros(M.size(0), dtype=N.complex128)
dofnos, rhs_contrib = point_source.calc_pointsource_contrib(
    V, source_coord, source_value)
b[dofnos] += rhs_contrib



# ml = smoothed_aggregation_solver(Asp,max_coarse=10)
# residuals = []
# x = ml.solve(b,tol=1e-10,accel='cg',residuals=residuals)

# residuals = residuals/residuals[0]
# print ml

