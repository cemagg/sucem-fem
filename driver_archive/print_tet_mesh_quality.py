"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
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
from NewCode.Utilities import Struct, partial
from NewCode import DifferentialForm, Waveforms, PostProc
from NewCode.DifferentialForm import Discretiser
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets
from NewCode.DiscretisedSystem import CurlCurlNewmark

mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh = mesh4

print 'Mesh elements: ', len(mesh.elements)

el = [N.linalg.norm(e.nodeCoords[0]-e.nodeCoords[1]) for e in mesh.edges]

print 'Edge length quality:'
print 'max: ', N.max(el),  '\tmin: ',  N.min(el), '\tavg: ', N.average(el)
rho_e = N.min(el)/N.max(el)
print 'rho_e: ', rho_e , '1/rho_e:', 1/rho_e

print 'Volume quality:'
vols = [elem.size for elem in mesh.elements]
print 'max: ', N.max(vols),  '\tmin: ',  N.min(vols), '\tavg: ', N.average(vols)
rho_v = N.min(vols)/N.max(vols)
print 'rho_v: ', rho_v , '1/rho_v:', 1/rho_v
