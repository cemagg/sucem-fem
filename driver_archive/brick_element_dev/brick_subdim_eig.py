"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
from itertools import izip
import os, sys
import pickle
import random
import numpy as N
import scipy
#
# Local Imports
#
import NewCode
from NewCode import PostProc
from NewCode.Meshes import BrickMesh
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import BrickSubDimMesh, DifferentialForm, SystemMatrix
from NewCode.DifferentialForm import BrickDiscretiser, BrickSubDimDiscretiser,\
     BrickSubDimDiscretiserEntities

import brick_cavity

a,b,l = 1., 1/4, 1/5
eps = 1e-10
edge1_x = close_to_point(0, eps)
edge2_x = close_to_point(a, eps)
edge3_y = close_to_point(b, eps)
edge4_y = close_to_point(0, eps)
z_port = close_to_point(0, eps)


mesh = BrickMesh.Mesh(brick_cavity.make_rect_cavity_brick_listmesh(
    a,b,l, [1/10, 1/10, 1/10]))

print 'Mesh elements: ', len(mesh.elements)

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

fs = lambda face: N.all(N.abs(face.nodeCoords[:,2]) < 0.000000001)
subFree =  lambda edge: not N.all([
    edge1_x(x) or edge2_x(x) or edge3_y(y) or edge4_y(y)
    for (x,y,z) in edge.nodeCoords])
direchFree = lambda ent: subFree(ent) and fs(ent)

subsurf = BrickSubDimMesh.SubSurface(mesh, fs)
disc = BrickDiscretiser.setup_PformDiscretiser(mesh, 1)
disc_direc = BrickDiscretiser.setup_PformDiscretiser(mesh,1,freeFun=direchFree)
subGeomEntities = {'edge': BrickSubDimDiscretiserEntities.SubDimDiscretiserEntityList(
    BrickSubDimDiscretiserEntities.Edge(mesh=subsurf, freefun=subFree,
                                        attrs=subsurf.edges.list_repr())),}
subDisc = BrickSubDimDiscretiser.PformSubDimDiscretiser(
    1, subsurf, subGeomEntities, BrickDiscretiser.Permuter, disc)

M = SystemMatrix.self_projection_matrix(subDisc)
S = SystemMatrix.self_projection_matrix(subDisc.D())
w,v = scipy.linalg.eig(S.todense(), M.todense())

rw = N.real(w)[N.abs(w) > 0.1]
rv = v.T[N.abs(w) > 0.1]

sort_ind = N.argsort(rw)

sw = rw[sort_ind[0:10]]
sv = rv[sort_ind[0:10]]

direch_dofs = disc_direc.newDOFs()
direch_dofs.dofArray[:] = sv[0]

modes = sorted([(m*N.pi/a)**2 + (n*N.pi/b)**2 for m in range(1,7) for n in range(0,7)])

cd = N.array([1,0,0], N.float64)
test_crossection = PostProc.MakeLine(0.01*a*cd, cd*a*0.99, 51)
reconst_pts = PostProc.LocatePoints(mesh, test_crossection)

r_field = direch_dofs.reconstruct(*reconst_pts)
plot(r_field[:,1])



