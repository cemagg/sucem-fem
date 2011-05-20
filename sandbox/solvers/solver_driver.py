__author__ = "Evan Lezar"
__date__ = "20 May 2011"

"""
A simple driver to test various LA solvers.
"""

from dolfin import dot, cross, curl, inner, dx, ds
import dolfin

import numpy as np
import os
import sys


sys.path.insert(0, '../../')
from FenicsCode.Consts import eps0, mu0, c0, Z0
from FenicsCode.Utilities.MeshIO import femmesh_2_dolfin_mesh
from FenicsCode.Utilities.Converters import dolfin_ublassparse_to_scipy_csr
from FenicsCode.Sources import point_source
from FenicsCode.Utilities.MatrixIO import ( load_scipy_matrix_from_mat, save_scipy_matrix_as_mat)
del sys.path[0]

def generate_matrices ( mesh_file ):
    """
    Generate the matrices for a 1\lambda spherical domain with an infintesimal dipole at the origin.
    
    @param mesh_file: the full path and filename to the spherical mesh to be used (this is a femmesh file).
    @todo: change the mesh type to a gmsh mesh.
    """
    problem_data = {'l': 2.9979245800000002e-4,
                    'I': 1.0,
                    'f': 1.0e9,
                    'source_coord': np.array([0,0,0.]),
                    'order': 1,
                    }
    
    lam = c0/problem_data['f']
    
    # load the mesh and scale coordinates accordingly
    mesh = femmesh_2_dolfin_mesh(mesh_file)
    mesh.init()
    mesh.coordinates()[:] *= lam
    
    
    eps_r = 1;
    mu_r = 1;
    freq = problem_data['f'];
    source_coord = problem_data['source_coord']
    source_point = dolfin.Point(*source_coord)
    source_value = np.array([0,0,1.])*problem_data['I']*problem_data['l']
    
    V = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", problem_data['order'])
    print 'DOFs: ', V.dim()

    # Define basis and bilinear form
    k_0 = 2*np.pi*freq/c0

    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)

    m = eps_r*inner(v, u)*dx                # Mass form
    s = (1/mu_r)*dot(curl(v), curl(u))*dx   # Stiffness form

    n = V.cell().n
    s_0 = inner(cross(n, v), cross(n, u))*ds

    def boundary(x, on_boundary):
        return on_boundary

    # Assemble forms
    print 'assembling forms'
    M = dolfin.uBLASSparseMatrix()
    S = dolfin.uBLASSparseMatrix()
    S_0 = dolfin.uBLASSparseMatrix()
    dolfin.assemble(m, tensor=M, mesh=mesh)
    dolfin.assemble(s, tensor=S, mesh=mesh)
    dolfin.assemble(s_0, tensor=S_0, mesh=mesh)

    # Set up RHS
    b = np.zeros(M.size(0), dtype=np.complex128)
    dofnos, rhs_contrib = point_source.calc_pointsource_contrib(
        V, source_coord, source_value)

    rhs_contrib = 1j*k_0*Z0*rhs_contrib
    b[dofnos] += rhs_contrib

    Msp = dolfin_ublassparse_to_scipy_csr(M)
    Ssp = dolfin_ublassparse_to_scipy_csr(S)
    S_0sp = dolfin_ublassparse_to_scipy_csr(S_0)
    A = Ssp - k_0**2*Msp + 1j*k_0*S_0sp

    return A, b

def save_matrices ( path, id, A, b ):
    """
    Save the system matrix and the excitation vector.
    The matrices will be saved in <path>/<id>/A.mat and <path>/<id>/b.mat
    
    @param path: the folder where the matrices are to be saved
    @param id: a unique problem identifier
    @param A: the system matrix
    @param b: the excitation vector
    """
    full_path = os.path.join(path, id)
    
    save_scipy_matrix_as_mat ( full_path, 'A', A )
    save_scipy_matrix_as_mat ( full_path, 'b', b )

def load_matrices ( path, id ):
    full_path = os.path.join ( path, id )
    A = load_scipy_matrix_from_mat ( full_path, 'A')
    b = load_scipy_matrix_from_mat ( full_path, 'b' )
    
    return A, b 


def generate_and_save ( problem_id ):
    mesh_file = '../../workspace/%s.femmesh' % problem_id
    data_path = 'matrices'
    A, b = generate_matrices ( mesh_file )
    save_matrices ( data_path, problem_id, A, b )
    return A, b
    
def load ( problem_id ):
    data_path = 'matrices'
    A, b = load_matrices ( data_path, problem_id )
    return A, b
    
def main ():
    problem_id = 'sphere-r1m-4'
    
    A_s, b_s = generate_and_save( problem_id )
    A, b = load ( problem_id )

    print b_s - b[:,0]
    
if __name__ == "__main__":
    main()