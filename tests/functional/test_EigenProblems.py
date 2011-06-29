__author__ = "Evan Lezar"
__date__ = "29 June 2011"

"""
This is a functional test for the solution of eigenproblems
"""

import dolfin
import numpy as np
import sys
import unittest

sys.path.insert(0, '../../')
from FenicsCode.ProblemConfigurations.EMVectorWaveEigen import EigenProblem
from FenicsCode.ProblemConfigurations.EMVectorWaveEigen import DefaultEigenSolver
from FenicsCode.Consts import c0
del sys.path[0]

def k_mnl ( abd, m, n, l=0, normalize = False):
    """
    Return the analytical cutoff wavenumber for a resonant cavity with dimensions specified in the tuple abd
    
    @param abd: a 3-tuple of the dimensions in meters of the cavitys
    @param m: the index of the mode in the x direction
    @param n: the index of the mode in the y direction
    @param l: the index of the mode in the z direction
    @param normalize: divide by the factor \pi
    """
    if len(abd) == 3:
        a, b, d = abd
    elif len(abd) == 2:
        a, b = abd
        d = 1;
        l = 0;
    
    root = np.sqrt ( (m/a)**2 + (n/b)**2 + (l/d)**2 )
    if normalize:
        return root
    else:
        return root*np.pi


class Test3D(unittest.TestCase):
    def test_cube_resonant_cavity(self):
        mesh = dolfin.UnitCube ( 4, 4, 4 )
        a = 1.0
        b = 1.0
        d = 1.0
        mesh.coordinates()[:,0] = a*mesh.coordinates()[:,0]
        mesh.coordinates()[:,1] = b*mesh.coordinates()[:,1]
        mesh.coordinates()[:,2] = d*mesh.coordinates()[:,2]
         
        # Use 3rd order basis functions 
        order = 3
        # Set up the eigen problem
        ep = EigenProblem()
        ep.set_mesh(mesh)
        ep.set_basis_order(order)
        ep.set_boundary_conditions(pec=True)
        ep.init_problem()
        
        # Set up eigen problem solver where sigma is the shift to use in the shift-invert process
        sigma = 1.1
        es = DefaultEigenSolver()
        es.set_eigenproblem(ep)
        es.set_sigma(sigma)
        
        # Solve the eigenproblem
        eigs_w, eigs_v = es.solve_problem(10)
        
        # Output the results
        res = np.array(sorted(eigs_w)[0:])
        res = np.sqrt(res)/np.pi
        
        steps = 3
        abd = (a, b, d)
        ids = []
        values = []
        for m in range(steps):
            for n in range(steps):
                for l in range(steps):
                    i = (m,n,l)
                    if i.count(0) < 2:
                        ids.append((m,n,l))
                        values.append(k_mnl ( abd, m, n, l, True ))
                    
                    # mode is both a TE and TM mode            
                    if i.count( 0 ) == 0:
                        ids.append((m,n,l))
                        values.append(k_mnl ( abd, m, n, l, True ))
        
        r = 0;
        errors = np.zeros_like(res)
        for i in np.argsort(values).tolist():
            if r < len(res):
                errors[r] = np.linalg.norm( res[r] - values[i])/np.linalg.norm( values[i] )
                r += 1
            else:
                break;
        
        np.testing.assert_array_almost_equal( errors, np.zeros_like(res), 4 )


class Test2D(unittest.TestCase):
    def test_TE_modes(self):
        mesh = dolfin.UnitSquare ( 5, 5 )
        a = 1.0
        b = 1.0
        mesh.coordinates()[:,0] = a*mesh.coordinates()[:,0]
        mesh.coordinates()[:,1] = b*mesh.coordinates()[:,1]
 
        # Use 3rd order basis functions 
        order = 3
        # Set up the eigen problem
        ep = EigenProblem()
        ep.set_mesh(mesh)
        ep.set_basis_order(order)
        ep.set_boundary_conditions(pec=True)
        ep.init_problem()
        
        # Set up eigen problem solver where sigma is the shift to use in the shift-invert process
        sigma = 1.5
        es = DefaultEigenSolver()
        es.set_eigenproblem(ep)
        es.set_sigma(sigma)

        # Solve the eigenproblem
        eigs_w, eigs_v = es.solve_problem(10)

        # Output the results
        res = np.array(sorted(eigs_w)[0:])
        res = np.sqrt(res)/np.pi
        
        steps = 5
        abd = (a, b)
        ids = []
        values = []
        for m in range(steps):
            for n in range(steps):
                l = 0
                i = (m,n)
                if i.count(0) < 2:
                    ids.append((m,n,l))
                    values.append(k_mnl ( abd, m, n, l, True ))

        r = 0;
        errors = np.zeros_like(res)
        for i in np.argsort(values).tolist():
            if r < len(res):
                errors[r] = np.linalg.norm( res[r] - values[i])/np.linalg.norm( values[i] )
                
                r += 1
            else:
                break;

        np.testing.assert_array_almost_equal( errors, np.zeros_like(res), 4 )
    
    def test_TM_modes(self):
        mesh = dolfin.UnitSquare ( 5, 5 )
        a = 1.0
        b = 1.0
        mesh.coordinates()[:,0] = a*mesh.coordinates()[:,0]
        mesh.coordinates()[:,1] = b*mesh.coordinates()[:,1]
        
        order = 4
        # Set up the eigen problem
        ep = EigenProblem()
        
        ep.set_mesh(mesh)
        ep.set_basis_order(order)
        ep.init_problem()
        
        # Set up eigen problem solver where sigma is the shift to use in the shift-invert process
        sigma = 1.5
        es = DefaultEigenSolver()
        es.set_eigenproblem(ep)
        es.set_sigma(sigma)
    
        # Solve the eigenproblem
        eigs_w, eigs_v = es.solve_problem(10)
    
        # Output the results
        res = np.array(sorted(eigs_w)[0:])
        res = np.sqrt(res)/np.pi
        
        steps = 5
        abd = (a, b)
        ids = []
        values = []
        for m in range(steps):
            for n in range(steps):
                i = (m,n)
                if i.count(0) == 0:
                    ids.append(i)
                    values.append(k_mnl ( abd, m, n, 0, True ))

        r = 0;
        errors = np.zeros_like(res)
        for i in np.argsort(values).tolist():
            if r < len(res):
                errors[r] = np.linalg.norm( res[r] - values[i])/np.linalg.norm( values[i] )
        
                r += 1
            else:
                break;
        
        np.testing.assert_array_almost_equal( errors, np.zeros_like(res), 4 )


if __name__ == "__main__":
    unittest.main()