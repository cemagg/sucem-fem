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
mesh=mesh4

#Discretiser.BasePformDiscretiser.defaultIntegrationOrder = 2

print 'Mesh elements: ', len(mesh.elements)

dt = 1/30.
cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

disc1f = Discretiser.setup_PformDiscretiser(mesh, 1, order=2, freeFun=free)
disc1f_curl = disc1f.D()
disc2f = Discretiser.setup_PformDiscretiser(mesh, 2, mixed=False,freeFun=free)
a,b,d,E0 = 0.5, 0.75, 1.0, 1.0
matchfield = lambda xyz: E0*N.array([
    0., 0., N.sin(N.pi*((xyz[1] - a/2)/a))*N.sin(N.pi*(xyz[0] - d/2)/d)])
curl_matchfield = lambda xyz: E0*N.array(
    [(N.pi*N.cos(N.pi*(xyz[1] - a/2)/a)*N.sin(N.pi*(xyz[0] - d/2)/d))/a,
     -(N.pi*N.sin(N.pi*(xyz[1] - a/2)/a)*N.cos(N.pi*(xyz[0] - d/2)/d))/d,
     0.])
E_dofs = disc1f.newDOFs()
print "Matching E-field"
E_dofs.matchFunction(matchfield)
print "Done"
E_dofs_curl = disc1f_curl.newDOFs()
E_dofs_curl.dofArray[:] = E_dofs.dofArray
match_curl_E_dofs = disc2f.newDOFs()
print "Matching curl(E-field)"
match_curl_E_dofs.matchFunction(curl_matchfield)
print "Done"
print "Calculating curl of matched E-field"
calc_curl_E_dofs = E_dofs.D(disc2f)
print "Calculating E-match error"
E_err = E_dofs.matchErrRMS(matchfield)
print "Calculating curl(E)-match error"
match_curl_err = match_curl_E_dofs.matchErrRMS(curl_matchfield)
print "Calculating calculated curl(E)-match error"
calc_curl_err = calc_curl_E_dofs.matchErrRMS(curl_matchfield)
print "Calculating curl(E)-match error of disc1f.D() version of E"
E_D_err = E_dofs_curl.matchErrRMS(curl_matchfield)

print "E_err"
print 100*E_err
print "match_curl_err"
print 100*match_curl_err
print "calc_curl_err"
print 100*calc_curl_err
print "E_D_err"
print 100*E_D_err
