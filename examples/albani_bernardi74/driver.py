from __future__ import division
__author = "Neilen Marais, Evan Lezar"
__date__ = "14 April 2011"

"""
This file has been adapted from a file eig_driver.py developed by N Marais

Problem considered:
    
    
Reference:
    Albani, Bernardi
"""
import sys
import numpy as N
import os
import dolfin as dol

sys.path.insert(0, '../../')
from FenicsCode.ProblemConfigurations.EMVectorWaveEigen import EigenProblem
from FenicsCode.ProblemConfigurations.EMVectorWaveEigen import DefaultEigenSolver
from FenicsCode.Consts import c0
del sys.path[0]

# Load the mesh and the material region markers
mesh_file = 'mesh/albani_bernardi74_fig2VII.xml'
materials_mesh_file = "%s_physical_region%s" % (os.path.splitext(mesh_file))
mesh = dol.Mesh(mesh_file)
material_mesh_func = dol.MeshFunction('uint', mesh, materials_mesh_file)
# Define the dielectric properties of the regions in the mesh
materials = {1000:dict(eps_r=16),
             1001:dict(eps_r=1)}

# Use 3rd order basis functions 
order = 3
# Set up the eigen problem
ep = EigenProblem()
ep.set_mesh(mesh)
ep.set_basis_order(order)
ep.set_material_regions(materials)
ep.set_region_meshfunction(material_mesh_func)
ep.init_problem()

# Set up eigen problem solver where sigma is the shift to use in the shift-invert process
sigma = 1.5
es = DefaultEigenSolver()
es.set_eigenproblem(ep)
es.set_sigma(sigma)

# Solve the eigenproblem
eigs_w, eigs_v = es.solve_problem()

# Output the results
res = N.array(sorted(eigs_w)[0:10])
print N.sqrt(res)
print c0*N.sqrt(res)/2/N.pi/1e6

print '\nreference:'
print "[ 2.40690704  2.59854088  3.2846475   3.85115539  4.03598938  4.2010115"
print "  4.38370524  4.48188875  4.78715748  4.86703817]"
print "[ 114.84184235  123.9853547   156.72186945  183.75191648  192.57098388"
print "  200.44475862  209.16170779  213.84638197  228.41180674  232.22319022]"
print '\n#TODO: replace with proper reference results'