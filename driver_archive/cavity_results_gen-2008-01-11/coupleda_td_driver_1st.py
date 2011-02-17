"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
import pickle
import numpy as N
import scipy
from scipy import signal
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
from NewCode.DiscretisedSystem import CoupledFirstOrderSystem

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

class TestRun(Runners.TestRun):
    drv_fun = staticmethod(Waveforms.get_gausspulse(fc, bw=0.9))
    useLU = True
    SystemClass = CoupledFirstOrderSystem

    def setupLogging(self):
        from NewCode import PostProc
        a,b,c = 29,23,19
        self.logPos = [ 17.,  14.,  12.]
        test_elnos, test_el_coords = PostProc.LocatePoints(
            self.system.discs.E.mesh, [self.logPos])
        self.system.addReconstructedLogger('E', test_elnos, test_el_coords)

    def getResult(self):
        tmp = N.array(self.system.loggedReconstructed['E'][0].vals).T
        return tmp.reshape(3,tmp.shape[2])

    def make_driveDOFs(self):
        print "Getting source DOFs"
        weights, elPerm = self.system.dofs.E.calcProjPointfunRHS_with_elPerm(
            matchfun=lambda r: N.array([1,1,1.]), r0=N.array([7.,2.,4.]))
        drive_dofnos = elPerm[1]
        weights = weights[elPerm[0]]
        return Struct(weights=weights, drive_dofnos=drive_dofnos)

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
#orders = (1,2,3,4)
#orders = (4,)
#orders = (2,3)
tms = (2500, 5000,10000,20000,40000,)#80000,160000)
# Timesteps smaller than ~0.0385 should be stable up to order 4 with
# rect-white-3 mesh
common_dts = (0.125, 0.0625, 0.03125, 0.015625)
# Largest stability timestep and then the power two divisors 0.5 smaller than
# stability step
initial_dts = {1: (3.18, 0.5, 0.25),
               2: (1.42, 0.5, 0.25),
               3: (0.7, 0.5, 0.25),
               4: (0.16)}
resfile = '+coupleda_results'

orders = (1,)                           # Choose desired orders here
for order in orders:
    times = dict((dt, tms) for dt in initial_dts[order] +
                 common_dts)
    E_order = H_order = (order, True)
    B_order = D_order = (order-1, False) if order > 1 else (order, True)

    runny = TestRun(mesh, BCs=Struct(E=cb, H=free, B=cb, D=free),
                    disc_orders=dict(E=E_order, B=B_order, H=H_order, D=D_order))
    for (dt, time), res in runny.yieldResults(times):
        pickle.dump({order:{(dt,time): res}}, file(
            resfile + '_o_%s_dt_%s_t_%s' % (
            str(order), str(dt), str(time)) + '.pickle', 'w'))



cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

