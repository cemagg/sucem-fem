__author__ = "Evan Lezar"
__date__ = "27 June 2011"

import sys
import numpy as N
import os
import dolfin as dol

sys.path.insert(0, '../../')
from FenicsCode.ProblemConfigurations.EMVectorWaveEigenproblem import EigenProblem, DefaultEigenSolver
from FenicsCode.Consts import c0
del sys.path[0]

script_path = os.path.dirname(__file__)
# Load the mesh and the material region markers
mesh = dol.UnitCube ( 4, 4, 4 )
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
    a, b, d = abd
    
    root = N.sqrt ( (m/a)**2 + (n/b)**2 + (l/d)**2 )
    if normalize:
        return root
    else:
        return root*N.pi
    
print '\nreference:'
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