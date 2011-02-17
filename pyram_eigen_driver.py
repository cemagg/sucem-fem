from __future__ import division
"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
import os
import pickle
from numpy import *
import numpy as N
import scipy
#
# Local Imports
#
import NewCode
from NewCode.Meshes import PyramMesh, BrickMeshGen
from NewCode import DifferentialForm, Utilities
from NewCode.Utilities import close_to_point
from NewCode.DifferentialForm import PyramDiscretiser
from NewCode.DifferentialForm import allfree
from NewCode.tests.PyramMeshes import TwelvePyram

#h0 = 1.0001/5.
a,b,c = 29,23,19
h = a/4.
#a,b,c = 1.,0.75,0.5
g_eps = 1e-10                           # Geometrical tollerance

zero_p = close_to_point(0, g_eps)
a_p, b_p, c_p = (close_to_point(xx, g_eps) for xx in (a,b,c))
def freeE(ent):
    x,y,z = ent.nodeCoords.T
    return not (N.all(a_p(x)) or N.all(zero_p(x)) or 
                N.all(b_p(y)) or N.all(zero_p(y)) or
                N.all(c_p(z)) or N.all(zero_p(z)))
mesh = PyramMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [h, h, h]))

#mesh = PyramMesh.Mesh(TwelvePyram.listmesh)
order = 2
mixed = True

print 'Mesh elements: ', len(mesh.elements)

class PyramVectorWaveEigen(NewCode.DiscretisedSystem.VectorWaveEigen):
    DiscretiserModule = PyramDiscretiser
    def __init__(self, mesh, order, mixed=True, BC=allfree,
                 volBC=None, btype=None):
        self.order = order
        self.mesh = mesh
        self.disc = self.DiscretiserModule.setup_PformDiscretiser(
            mesh, self.p, order, mixed, BC, volBC, btype)
        self.disc.setIntegrationRule((3,4))
        self.disc.D().setIntegrationRule((3,4))
        self.massMatrix = self.disc.matrix.mass
        self.stiffnessMatrix = self.disc.matrix.stiffness


# eigsys = PyramVectorWaveEigen(mesh, order, mixed=mixed, BC=freeE,
#                               btype='graglia99')

eigsys = PyramVectorWaveEigen(mesh, order, mixed=mixed, BC=freeE,
                              btype='coulomb97')

print 'Getting S'
S = eigsys.stiffnessMatrix()
print 'Getting M'
M = eigsys.massMatrix()

from scipy.sparse.linalg.eigen.arpack import speigs
sigma = 0.01
print "(nodofs, nnz, sparsity %)", M.shape[0], M.nnz, M.nnz/M.shape[0]**2*100
# print 'Sparse LU decomposition'
# sigma_solve = scipy.sparse.linalg.dsolve.factorized(S - sigma*M)

print 'Solving Eigenproblem'
# w,v = speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, M.shape[0], sigma, 51, ncv=91)
w,v = scipy.linalg.eig(S.todense(), M.todense())
from AnalyticResults import PEC_cavity, accoustic_eigs, err_percentage

ares = PEC_cavity['rect-white']
#ares = PEC_cavity['rect1x0.75x0.5']
res =  N.array(sorted(N.abs(w[w > 0.0000001]))[0:10])
err = err_percentage(ares, res)
RMS_err = Utilities.RMS(err[0:4])
#edges = eigsys.disc.geomEntities['edge']
#faces = eigsys.disc.geomEntities['face']

el = eigsys.disc.elements[0]
el_D = eigsys.disc.D().elements[0]

print RMS_err
print err
