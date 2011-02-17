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
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import SubDimMesh, DifferentialForm, SystemMatrix, PostProc
from NewCode.DifferentialForm import Discretiser, SubDimDiscretiser,SubDimDiscretiserEntities
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets
from NewCode import Feeds

a,b,l = 1., 1/4, 1/5

mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh=mesh4

print 'Mesh elements: ', len(mesh.elements)

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

fs = lambda face: N.all(N.abs(face.nodeCoords[:,2]) < 0.000000001)

eps = 1e-10
edge1_x = close_to_point(0, eps)
edge2_x = close_to_point(a, eps)
edge3_y = close_to_point(b, eps)
edge4_y = close_to_point(0, eps)
z_port = close_to_point(0, eps)

subFree =  lambda edge: not N.all([
    edge1_x(x) or edge2_x(x) or edge3_y(y) or edge4_y(y)
    for (x,y,z) in edge.nodeCoords])
#subFree = free

direchFree = lambda ent: subFree(ent) and fs(ent)

bf_order = 1

subsurf = SubDimMesh.SubSurface(mesh, fs)
disc = Discretiser.setup_PformDiscretiser(mesh, 1, bf_order)
disc_direc = Discretiser.setup_PformDiscretiser(mesh,1,bf_order, freeFun=direchFree)
subGeomEntities = {'edge': SubDimDiscretiserEntities.SubDimDiscretiserEntityList(
    SubDimDiscretiserEntities.Edge(mesh=subsurf, freefun=subFree,
                                   attrs=subsurf.edges.list_repr()))}
if bf_order >= 2:
    subGeomEntities['face'] =  SubDimDiscretiserEntities.SubDimDiscretiserEntityList(
        SubDimDiscretiserEntities.Face(mesh=subsurf, freefun=subFree,
                                       attrs=subsurf.elements.list_repr()))
subDisc = SubDimDiscretiser.PformSubDimDiscretiser(
    1, subsurf, subGeomEntities, Discretiser.Permuter, disc)

M = SystemMatrix.self_projection_matrix(subDisc)
S = SystemMatrix.self_projection_matrix(subDisc.D())
w,v = scipy.linalg.eig(S.todense(), M.todense())

rw = N.real(w)[N.abs(w) > 0.1]
rv = v.T[N.abs(w) > 0.1]

sort_ind = N.argsort(rw)

sw = rw[sort_ind[0:10]]
sv = rv[sort_ind[0:10]]

direch_dofs = disc_direc.newDOFs()
direch_dofs.dofArray[0:len(sv[0])] = sv[0]


modes = sorted([(m*N.pi/a)**2 + (n*N.pi/b)**2
                for m in range(0,7) for n in range(0,7)
                if (m > 0) or (n > 0)])[0:10]

cd = N.array([1,0,0], N.float64)
test_crossection = PostProc.MakeLine(0.01*a*cd, cd*a*0.99, 51)
reconst_pts = PostProc.LocatePoints(mesh, test_crossection)

r_field = direch_dofs.reconstruct(*reconst_pts)

#plot(r_field[:,1])

wgf = Feeds.WaveguideEigenMatcher(disc)
wgf_w, wgf_v = wgf.solve(10)
