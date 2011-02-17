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
from scipy import linalg
#
# Local Imports
#
import NewCode
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Utilities import Struct, partial
from NewCode import DifferentialForm, Waveforms, PostProc
from NewCode.DifferentialForm import Discretiser
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets
from NewCode.DiscretisedSystem import CurlCurlNewmark, CoupledFirstOrderSystem


mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh=mesh4

print 'Mesh elements: ', len(mesh.elements)

Discretiser.BasePformDiscretiser.defaultIntegrationOrder = 8

weights = 1.0
dt = 0.01
cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
# mixed 4th
orders = {'E':(4,True), 'H':(4,True), 'D':(3,False), 'B':(3,False)}
# mixed 3rd
#orders = {'E':(3,True), 'H':(3,True), 'D':(2,False), 'B':(2,False)}
# LT/QN E, H; Linear D, B
#orders = {'E':(2,True), 'H':(2,True), 'D':(1,False), 'B':(1,False)}
# All lowest order
#orders = {'E':(1,True), 'H':(1,True), 'D':(1,True), 'B':(1,True)}
coupledSystem = CoupledFirstOrderSystem(
    mesh, BCs=Struct(E=cb, H=free, D=free, B=cb),
    disc_orders=orders)

totalDOFs_D = coupledSystem.discs.D.totalDOFs
totalDOFs_E = coupledSystem.discs.E.totalDOFs

E,D,H,B = [coupledSystem.discs[dn] for dn in ('E', 'D', 'H', 'B')]

a,b,c = 29.,23.,19.
matchfun = lambda r: N.array([0,0,N.sin(N.pi*r[1]/b)*N.sin(N.pi*r[0]/a)],
                             N.float64)

match_line = PostProc.MakeLine(N.array([0,b/2,c/2]), N.array([a,b/2,c/2]), 101)
dr = 0.25
match_surf_x = int(a/dr)
match_surf_y = int(b/dr)
match_surf_line = PostProc.MakeLine(N.array([0,0,c/2]), N.array([a,0,c/2]), match_surf_x)
match_surf = N.array([match_surf_line + [[0,y,0]] for y in N.arange(match_surf_y)*dr],
                     N.float64).reshape(match_surf_x*match_surf_y,3)

elnos, local_coords = PostProc.LocatePoints(E.mesh, match_line)
elnos_surf, local_coords_surf = PostProc.LocatePoints(E.mesh, match_surf)
Edofs = E.newDOFs()
eig_modes = pickle.load(file('+coupleda_4th_eigmodes.pickle'))
Edofs.dofArray[:] = eig_modes[0.030392653822920439]
recr_field = Edofs.reconstruct(elnos, local_coords)

recr_field_surface_EBHD_good = Edofs.reconstruct(elnos_surf, local_coords_surf)
Edofs.dofArray[:] = eig_modes[0.030297962050100127]
recr_field_surface_EBHD_bad = Edofs.reconstruct(elnos_surf, local_coords_surf)
EB_eig_modes = pickle.load(file('+coupledb_4th_eigmodes.pickle'))
Edofs.dofArray[:] = EB_eig_modes[0.030392655151965684]
recr_field_surface_EB = Edofs.reconstruct(elnos_surf, local_coords_surf)
analytic_field_surface = N.array([matchfun(r) for r in match_surf], N.float64)
subplot(221)
pcolor(recr_field_surface_EBHD_good[:,2].reshape(match_surf_y,match_surf_x),
       cmap=cm.gray)
axis('tight')
xticks([])
yticks([])
title('EBHD mixed 4th, k=0.1743349')
subplot(222)
pcolor(recr_field_surface_EBHD_bad[:,2].reshape(match_surf_y,match_surf_x),
       cmap=cm.gray)
axis('tight')
xticks([])
yticks([])
title('EBHD mixed 4th, k=0.1740631')
subplot(223)
pcolor(recr_field_surface_EB[:,2].reshape(match_surf_y,match_surf_x),
       cmap=cm.gray)
axis('tight')
xticks([])
yticks([])
title('EB mixed 4th, k=0.1743349')
subplot(224)
pcolor(analytic_field_surface[:,2].reshape(match_surf_y,match_surf_x),
       cmap=cm.gray)
axis('tight')
xticks([])
yticks([])
title('Analytical, k=0.1743349')

