from __future__ import division

import sys
sys.path.append('../../')
import numpy as N
import os
import dolfin as dol
from dolfin import dot, cross, curl, inner, dx, ds
import scipy.sparse

from FenicsCode.Consts import eps0, mu0, c0, Z0
from FenicsCode.Utilities.Converters import dolfin_ublassparse_to_scipy_csr
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
    # There seems to be an issue with the intersect operator if the
    # mesh coordinates are changed after calling it for the first
    # time. Since we called it to find the centroid, we should init a
    # new mesh
    mesh_coords = mesh.coordinates().copy()
    mesh = dol.UnitCube(*parameters['domain_subdivisions'])
    mesh.coordinates()[:] = mesh_coords
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

    rhs_contrib = 1j*k_0*Z0*rhs_contrib
    b[dofnos] += rhs_contrib


    import pyamg 
    
    print M.size(0)
    Msp = dolfin_ublassparse_to_scipy_csr(M)
    Ssp = dolfin_ublassparse_to_scipy_csr(S)
    S_0sp = dolfin_ublassparse_to_scipy_csr(S_0, dtype=N.complex128, imagify=True)

    A = Ssp - k_0**2*Msp + k_0*S_0sp 
    
    import scipy.sparse.linalg
    A_lu = scipy.sparse.linalg.factorized(A.T)
    x = A_lu(b)

    #print ml

    workspace['V'] = V
    workspace['u'] = u
    workspace['x'] = x
    workspace['A'] = A
    workspace['b'] = b

def get_E_field(workspace, field_pts):
    dol.set_log_active(False)
    x = workspace['x']
    V = workspace['V']
    mesh = V.mesh()
    u_re = dol.Function(V)
    u_im = dol.Function(V)
    u_re.vector()[:] = N.real(x)
    u_im.vector()[:] = N.imag(x)
    E_field = N.zeros((len(field_pts), 3), dtype=N.complex128)
    for i, fp in enumerate(field_pts):
        try: E_field[i,:] = u_re(fp) + 1j*u_im(fp)
        except RuntimeError: E_field[i,:] = N.nan + 1j*N.nan
    return E_field
