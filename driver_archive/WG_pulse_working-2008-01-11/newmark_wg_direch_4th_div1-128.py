"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
import pickle
import random
import numpy as N
import scipy
from scipy import integrate

#
# Local Imports
#
import NewCode
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import SubDimMesh, DifferentialForm, Waveforms, PostProc, Feeds
from NewCode.DifferentialForm import Discretiser, SubDimDiscretiser,SubDimDiscretiserEntities
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets
from NewCode.DiscretisedSystem import CurlCurlNewmarkDirechlet

from NewCode.Feeds import WaveguideEigenMatcher

mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh = mesh4

print 'Mesh elements: ', len(mesh.elements)


Discretiser.BasePformDiscretiser.defaultIntegrationOrder = 8

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

freeE = cb

from analytical_WG_driver import WGT
#WGT.calc_ts()
#WGT.calc_FFTs()

g_eps = 1e-10                           # Geometrical tollerance

drv_fun = WGT.discrete_drv_fn
z_measure = WGT.test_z
z_measure_p = close_to_point(z_measure, g_eps)
runtime = WGT.t_final
z_port = 0.
z_port_p = close_to_point(z_port, g_eps)
a=1. ; b=0.25                           # WG x/y dim
zero_p = close_to_point(0, g_eps)
a_p, b_p = (close_to_point(xx, g_eps) for xx in (a,b))
on_port = lambda ent: N.all(z_port_p(ent.nodeCoords[:,2]))
port_free = lambda ent: not N.all([ a_p(x) or zero_p(x) or b_p(y) or zero_p(y)
                                     for (x,y,z) in ent.nodeCoords])
on_measurement_port = lambda ent: N.all(z_measure_p(ent.nodeCoords[:,2]))
measurement_port_free = port_free
direch_free = lambda ent: on_port(ent) and port_free(ent)  


# cd = N.array([1,0,0], N.float64)
# test_crossection = PostProc.MakeLine(0.01*cd, cd*a*0.99, 51) + [0,0,0.3]
# newmarkSystem.addReconstructedLogger(
#     *PostProc.LocatePoints(newmarkSystem.disc.mesh, test_crossection))


class TestRun(object):
    drv_fun = staticmethod(drv_fun)
    useLU = True
    SystemClass = CurlCurlNewmarkDirechlet
    log_divisor = 1
    def __init__(self, mesh, **kwargs):
        self.system = self.SystemClass(mesh, **kwargs)

    def setupSource(self):
        self.system.drive_fun = self.drv_fun
        self.system.setDirechBCs(direch_free)
        self.sm = sm = Feeds.SurfaceFieldMatcher()
        sm.initSubdim(self.system.disc, on_port, port_free)
        self.E_wg = E_wg = Feeds.gen_E_TE01(a,b)
        self.direch_dofArray = sm.matchKnown(E_wg)
        self.system.direchSys.dofs.dofArray[:] = self.direch_dofArray

    def set_dt(self, dt):
        self.dt = dt
        self.resetHistory()
        if self.useLU: self.system.useLU = True
        print "Setting timestep"
        self.system.setTimestep(dt)
        print "Done"

    def get_stepsCompleted(self):
        return self.system.n

    def resetHistory(self):
        self.system.reset_history()
        self.setupLogging()
        self.system.log()                     # To add an entry for n=0
        
    def setupLogging(self):
        divisor = self.log_divisor
        from NewCode import PostProc
        self.sm_m = sm_m = Feeds.SurfaceFieldMatcher()
        sm_m.initSubdim(self.system.disc, on_measurement_port, measurement_port_free)
        self.measure_dofArray = sm_m.matchKnown(self.E_wg)
        self.tif = tif = Feeds.TangentialFuncProjSurfaceIntegral(
            self.system.disc, on_measurement_port, measurement_port_free)
        self.test_pts = test_pts = N.array([[a/2, b/2, z_measure]], N.float64)
        test_elnos, test_el_coords = PostProc.LocatePoints(self.system.disc.mesh, test_pts)
        self.test_elnos = test_elnos ; self.test_el_coords = test_el_coords
        self.system.addReconstructedLogger(test_elnos, test_el_coords, divisor=divisor)
        self.logged_dofnos = logged_dofnos = tif.superDOFMap
        self.system.addLogger(dofnos=logged_dofnos, divisor=divisor)
        
    def runSteps(self, n):
        self.system.step(n)

    def getResult(self):
        rsv = N.array(self.system.loggedReconstructed[0]['vals'])[:,0,:]
        point_reconstructed_ts = rsv[:,1]
        sm_m = self.sm_m
        measure_dofArray = self.measure_dofArray
        nf =  N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(measure_dofArray))
        ts_modeintg_n = N.array([
            N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(dof_vec))
            for dof_vec in self.system.loggedDOFs[tuple(self.logged_dofnos)].vals],
                                N.float64)/nf        
        return Struct(point_reconstructed_ts=point_reconstructed_ts,
                      nf=nf, ts_modeintg_n=ts_modeintg_n)

analytical_dt = WGT.dt                  # Should be 0.0625

#orders = (4,3,2,1)
orders = (4,)
dt_divisors = (128,64,32,16,8,4,2,1)
resses = {}
for order in orders:
    if not order in resses: resses[order]={}
    TRS = TestRun(mesh, order=order, BC=freeE, useQ=False)
    TRS.setupSource()
    for dt_div in dt_divisors:
        dt = analytical_dt/dt_div
        no_steps = int(N.ceil(runtime/dt))
        TRS.log_divisor = dt_div
        TRS.set_dt(dt)      
        TRS.runSteps(no_steps)
        resses[order][dt_div] = TRS.getResult()


pickle.dump(resses, file('+newmark_wg_direch_4th_div1-128.pickle', 'w'))

from analytical_WG_driver import n_f_min, n_f_max
#WGT.ts = res.ts_modeintg_n
#WGT.calc_FFTs()

#x = N.arange(n_f_min, n_f_max)*WGT.df
#plot(x, 20*N.log10(N.abs(WGT.transfer_fs[n_f_min:n_f_max])))


