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
from NewCode import DifferentialForm, Waveforms
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

Discretiser.BasePformDiscretiser.defaultIntegrationOrder = 8

from NewCode import Runners

fc = 0.0375
d_gaus = Waveforms.get_d_gaussian(fc=fc)

class TestRun(Runners.TestRun):
    #drv_fun = staticmethod(Waveforms.get_gausspulse(0.0375, bw=0.9))
    drv_fun = staticmethod(d_gaus.D)
    useLU = True
    SystemClass = CurlCurlNewmark

    def setupLogging(self):
        from NewCode import PostProc
        a,b,c = 29,23,19
        self.logPos = [ 17.,  14.,  12.]
        test_elnos, test_el_coords = PostProc.LocatePoints(
            self.system.disc.mesh, [self.logPos])
        self.system.addReconstructedLogger(test_elnos, test_el_coords)

    def getResult(self):
        tmp = N.array(self.system.loggedReconstructed[0].vals).T
        return tmp.reshape(3,tmp.shape[2])

        

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

#orders = (1,2,3,4)
#orders = (4,)
#orders = (2,3)
tms = (2500,5000,10000,20000,40000,80000)#,160000)
#times = {0.5:tms, 0.25:tms, 0.125:tms, 0.0625:tms, 0.03125:tms, 0.015625:tms,
#         0.0078125:tms }
times = {0.25:tms,}         
results = dict(info="Newmark beta using white97's geometry, tets, beta = 0.25001"
               +str(times))
resfile = '../working/+newmark_results_beta_0.25001'

orders = (4,)                           # Choose desired orders here
for order in orders:
    runny = TestRun(mesh, BC=cb, order=order)
    for (dt, time), res in runny.yieldResults(times):
        pickle.dump({order:{(dt,time): res}}, file(
            resfile + '_o_%s_dt_%s_t_%s' % (str(order), str(dt), str(time)) + '.pickle', 'w'))





