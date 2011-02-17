"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
import pickle
import numpy as N
import scipy
import random
from scipy import signal
#
# Local Imports
#
import NewCode
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Utilities import Struct, partial
from NewCode import DifferentialForm, Waveforms
from NewCode.DifferentialForm import Discretiser, DiscretiserDOFs
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets
from NewCode.DiscretisedSystem import CoupledFirstOrderSystemB

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
    drv_fun = staticmethod(d_gaus)
    useLU = True
    SystemClass = CoupledFirstOrderSystemB

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
orders = (1,2,3,4)
#orders = (4,)
#orders = (2,3)
tms = (2500,5000,10000,20000,40000,80000,160000)
# Timesteps smaller than ~0.0385 should be stable up to order 4 with
# rect-white-3 mesh
common_dts = (0.03125, 0.015625, 0.0078125)
# Largest stability timestep and then the power two divisors 0.5 smaller than
# stability step
initial_dts = {1: (1.2495, 0.5, 0.25, 0.125, 0.0625),
               2: (0.4288, 0.25, 0.125, 0.0625),
               3: (0.1895, 0.125, 0.0625),
               4: (0.0384,)}

orders = (1,)
results = dict(info="Leapfrog using white97's geometry, rect-white-3")
for order in orders:
    times = dict((dt, tms) for dt in initial_dts[order] + common_dts)
    E_order = (order, True)
    B_order = (order-1, False) if order > 1 else (order, True)
    runny = TestRun(mesh, BCs=Struct(E=cb, B=cb),
                    disc_orders=dict(E=E_order, B=B_order))
    results[order] = runny.getResults(times)

import pickle
pickle.dump(results, file('+coupledb_results_1st.pickle', 'w'))
