"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
import pickle
import numpy as N
import scipy
#
# Local Imports
#
import NewCode
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Utilities import Struct, partial
from NewCode import DifferentialForm
from NewCode.DifferentialForm import Discretiser
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets
from NewCode import MaximaUtils

mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh = mesh4

print 'Mesh elements: ', len(mesh.elements)

disc = Discretiser.setup_PformDiscretiser(mesh, 1, order=4)
disc_dofs = disc.newDOFs(N.complex128)
disc_dofs.matrix.dtype = disc_dofs.dtype
k = N.pi*2
matchfun = lambda r: N.array([N.e**(-1j*k*r[2]), 0, 0], N.complex128)
disc_dofs.matchFunction(matchfun)
r = 0.5
sphere_vol = 4/3.*N.pi*r**3
RMS_err = disc_dofs.matchErrRMS(matchfun)/sphere_vol
RMS_err_percent = RMS_err*100

# f = file('/tmp/exportmesh.mac', 'w')
# MaximaUtils.write_mesh(mesh, f)
# MaximaUtils.write_dofs(disc_dofs.dofArray, f)
# MaximaUtils.write_dofmap(disc, f)
# f.close()

### Results

"""
RMS err CT/LN Galerkin matching:

lam/4  : 8.9%
lam/8  : 6.1%
lam/16 : 3.2%
lam/32 : 1.6%

RMS err LT/QN Galerkin matching:

lam/4 : 1.6%?

RMS err (3,2), Galerkin:

lam/4 : 0.34%

RMS err (4,3), Galerkin:

lam/4 : 0.02%

"""
