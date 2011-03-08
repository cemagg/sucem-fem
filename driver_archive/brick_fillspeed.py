"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
import pickle
from numpy import *
import numpy as N
import scipy
#
# Local Imports
#
import NewCode
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode import DifferentialForm, Utilities
from NewCode.DifferentialForm import BrickDiscretiser
from NewCode.DifferentialForm import constrained_on_boundary, allfree

#h0 = 1/8.
a,b,c = 29,23,19
#a,b,c = 1.,0.75,0.5
#a,b,c = 1, 0.25, 5.1
h = a/10.
# mesh = BrickMesh.Mesh(
#     BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [a*h0, b*h0, c*h0]))

mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [h, h, h]))

order = 3

freeE = constrained_on_boundary

print 'Mesh elements: ', len(mesh.elements)

class BrickVectorWaveEigen(NewCode.DiscretisedSystem.VectorWaveEigen):
    DiscretiserModule = BrickDiscretiser
    def __init__(self, mesh, order, mixed=True, BC=constrained_on_boundary,
                 volBC=None):
        self.order = order
        self.mesh = mesh
        self.disc = self.DiscretiserModule.setup_PformDiscretiser(
            mesh, self.p, order, mixed, BC, volBC, 'cohen98')
        self.massMatrix = self.disc.matrix.mass
        self.stiffnessMatrix = self.disc.matrix.stiffness

class BrickAudioEigen(NewCode.DiscretisedSystem.AudioEigen):
    DiscretiserModule = BrickDiscretiser
    def __init__(self, mesh, order, mixed=True, BC=constrained_on_boundary,
                 volBC=None):
        self.order = order
        self.mesh = mesh
        self.disc = self.DiscretiserModule.setup_PformDiscretiser(
            mesh, self.p, order, mixed, BC, volBC, 'cohen98')
        self.massMatrix = self.disc.matrix.mass
        self.stiffnessMatrix = self.disc.matrix.stiffness

eigsys = BrickVectorWaveEigen(mesh, order, mixed=True, BC=freeE)

eigsys.disc.diagonalise()

# print 'Getting S'
# S = eigsys.stiffnessMatrix()
print 'Getting M'
M = eigsys.massMatrix()

# import cProfile
# cProfile.run('eigsys.massMatrix()', filename="+brick_filspeed.cprof", sort=1)
# import pstats
# stats = pstats.Stats("+brick_filspeed.cprof")
# stats.strip_dirs().sort_stats('cum', 'time').print_stats(20)

