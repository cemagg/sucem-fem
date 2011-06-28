__author__ = "Evan Lezar"
__date__ = "28 June 2011"

"""A simple 2D eigenproblem which calculates the TM modes of a square guide.

Note that this is done by modelling the Magnetic field and as such no dirichlet BCs are used.

Only natural boundary conditions."""

import sys
import numpy as N
import os
import dolfin as dol

sys.path.insert(0, '../../../')
from FenicsCode.ProblemConfigurations.EMVectorWaveEigen import EigenProblem
from FenicsCode.ProblemConfigurations.EMVectorWaveEigen import DefaultEigenSolver
from FenicsCode.Consts import c0
del sys.path[0]

script_path = os.path.dirname(__file__)
# Load the mesh and the material region markers
mesh = dol.UnitSquare ( 5, 5 )
a = 1.0
b = 1.0
mesh.coordinates()[:,0] = a*mesh.coordinates()[:,0]
mesh.coordinates()[:,1] = b*mesh.coordinates()[:,1]
 
# Use 3rd order basis functions 
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
res = N.array(sorted(eigs_w)[0:])
res = N.sqrt(res)/N.pi
print res



def k_mnl ( abd, m, n, l, normalize = False):
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
    
    root = N.sqrt ( (m/a)**2 + (n/b)**2 + (l/d)**2 )
    if normalize:
        return root
    else:
        return root*N.pi
    
print '\nreference:'
steps = 5
abd = (a, b)
ids = []
values = []
for m in range(steps):
    for n in range(steps):
        l = 0
        i = (m,n)
        if i.count(0) == 0:
            ids.append((m,n,l))
            values.append(k_mnl ( abd, m, n, l, True ))

r = 0;
errors = N.zeros_like(res)
print "mnl, analytical, calculated, relative error"
for i in N.argsort(values).tolist():
    if r < len(res):
        errors[r] = N.linalg.norm( res[r] - values[i])/N.linalg.norm( values[i] )
        print "%d%d%d, " % (ids[i]), "%9.3f, %10.3f, %.2e" % ( values[i], res[r], errors[r] )
        
        r += 1
    else:
        break;

N.testing.assert_array_almost_equal( errors, N.zeros_like(res), 4 )