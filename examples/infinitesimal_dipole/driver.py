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
# parameters dictionary, should be rebound by user of module
parameters = dict(f=None, l=None, I=None, source_coord=None)


def run(parameters, workspace):
    # Define material and frequency
    eps_r = 1;
    mu_r = 1;
    freq = parameters['f'];
    source_coord = parameters['source_coord']
    source_point = dol.Point(*source_coord)
    source_value = N.array([0,0,1.])*parameters['I']*parameters['l']

    # Define mesh
    mesh = dol.UnitCube(*parameters['domain_subdivisions'])
    # Transform mesh to correct dimensions
    mesh.coordinates()[:] *= parameters['domain_size']
    mesh.coordinates()[:] -= parameters['domain_size']/2
    ## Translate mesh slightly so that source coordinate lies at centroid of an element
    source_elnos = mesh.all_intersected_entities(source_point)
    closest_elno = source_elnos[(N.argmin([source_point.distance(dol.Cell(mesh, i).midpoint())
                                      for i in source_elnos]))]
    centre_pt = dol.Cell(mesh, closest_elno).midpoint()
    centre_coord = N.array([centre_pt.x(), centre_pt.y(), centre_pt.z()])
    mesh.coordinates()[:] -= centre_coord
    ##

    # Define function space
    order = parameters['order']
    V = dol.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", order)

    # Define basis and bilinear form
    k_0 = 2*N.pi*freq/c0

    u = dol.TrialFunction(V)
    v = dol.TestFunction(V)

    m = eps_r*inner(v, u)*dx                # Mass form
    s = (1/mu_r)*dot(curl(v), curl(u))*dx   # Stiffness form

    n = V.cell().n
    s_0 = inner(cross(n, v), cross(n, u))*ds

    def boundary(x, on_boundary):
        return on_boundary

    # Assemble forms
    M = dol.uBLASSparseMatrix()
    S = dol.uBLASSparseMatrix()
    S_0 = dol.uBLASSparseMatrix()
    dol.assemble(m, tensor=M, mesh=mesh)
    dol.assemble(s, tensor=S, mesh=mesh)
    dol.assemble(s_0, tensor=S_0, mesh=mesh)

    # Set up RHS
    b = N.zeros(M.size(0), dtype=N.complex128)
    dofnos, rhs_contrib = point_source.calc_pointsource_contrib(
        V, source_coord, source_value)
    b[dofnos] += rhs_contrib


    from scipy.sparse import csr_matrix
    import pyamg 
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

    residuals = []
    ml = pyamg.smoothed_aggregation_solver(A,max_coarse=10)
#    x = ml.solve(b,tol=1e-10,residuals=residuals)
#    import pdb; pdb.set_trace()
    x = ml.solve(b,tol=1e-10,accel='cg',residuals=residuals)
    # x, ml = pyamg.solve(A,b,verb=False, return_solver=True)
    import scipy.sparse.linalg
    A_lu = scipy.sparse.linalg.factorized(A)
    x_lu = A_lu(b)

    residuals = residuals
    #print ml

    workspace['V'] = V
    workspace['x'] = x
    workspace['x_lu'] = x_lu
    workspace['A'] = A
    workspace['ml'] = ml
    workspace['b'] = b
    workspace['residuals'] = N.array(residuals)
