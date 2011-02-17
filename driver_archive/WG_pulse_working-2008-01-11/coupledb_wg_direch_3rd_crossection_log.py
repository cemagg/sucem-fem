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
from NewCode import SubDimMesh, DifferentialForm, SystemMatrix, PostProc, Runners, Feeds
from NewCode.DifferentialForm import Discretiser, SubDimDiscretiser,SubDimDiscretiserEntities
from NewCode.DifferentialForm import DiscretiserDOFs
from NewCode.DiscretisedSystem import CoupledFirstOrderSystemBDirechlet
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets

from NewCode.Feeds import WaveguideEigenMatcher

mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh=mesh4

print 'Mesh elements: ', len(mesh.elements)

Discretiser.BasePformDiscretiser.defaultIntegrationOrder = 8

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
constrained = DifferentialForm.allconstrained
 
from analytical_WG_driver import WGT

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
freeE = cb
freeB = lambda ent: freeE(ent) or on_port(ent)
direch_free = lambda ent: on_port(ent) and port_free(ent)  

class TestRun(Runners.TestRun):
    drv_fun = staticmethod(drv_fun)
    useLU = True
    SystemClass = CoupledFirstOrderSystemBDirechlet

    def __init__(self, mesh, **kwargs):
        self.system = self.SystemClass(mesh, **kwargs)
        
    def setupSource(self):
        self.system.drive_fun = self.drv_fun
        self.sm = sm = Feeds.SurfaceFieldMatcher()
        sm.initSubdim(self.system.discs.E, on_port, port_free)
        self.E_wg = E_wg = Feeds.gen_E_TE01(a,b)
        self.direch_dofArray = sm.matchKnown(E_wg)
        self.system.setDirechBCs(Struct(E=direch_free, B=constrained))
        self.system.direchSys.dofs.E.dofArray[:] = self.direch_dofArray
        
    def setupLogging(self):
        divisor = self.log_divisor
        from NewCode import PostProc
        self.sm_m = sm_m = Feeds.SurfaceFieldMatcher()
        sm_m.initSubdim(self.system.discs.E, on_measurement_port, measurement_port_free)
        self.measure_dofArray = sm_m.matchKnown(self.E_wg)
        self.tif = tif = Feeds.TangentialFuncProjSurfaceIntegral(
            self.system.discs.E, on_measurement_port, measurement_port_free)
        self.test_pts = test_pts = N.array([[a/2, b/2, z_measure]], N.float64)
        test_elnos, test_el_coords = PostProc.LocatePoints(self.system.discs.E.mesh, test_pts)
        self.test_elnos = test_elnos ; self.test_el_coords = test_el_coords
        self.system.addReconstructedLogger('E', test_elnos, test_el_coords, divisor=divisor)
        self.logged_dofnos = logged_dofnos = tif.superDOFMap
        self.system.addLogger('E', dofnos=logged_dofnos, divisor=divisor)
        ## for cross-ection logging
        cd = N.array([a,0,0], N.float64)
        self.crossec_pts = crossec_pts = PostProc.MakeLine(
            0*cd, cd, 101) + [0,b/2.,z_measure]
        crossec_elnos, crossec_coords = PostProc.LocatePoints(
            self.system.discs.E.mesh, crossec_pts)
        self.system.addReconstructedLogger('E', crossec_elnos, crossec_coords, divisor=divisor)

    def getResult(self):
        rsv = N.array(self.system.loggedReconstructed.E[0]['vals'])[:,0,:]
        point_reconstructed_ts = rsv[:,1]
        sm_m = self.sm_m
        measure_dofArray = self.measure_dofArray
        nf =  N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(measure_dofArray))
        ts_modeintg_n = N.array([
            N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(dof_vec))
            for dof_vec in self.system.loggedDOFs.E[tuple(self.logged_dofnos)].vals],
                                N.float64)/nf        
        test_crossec_field = N.array(self.system.loggedReconstructed.E[1]['vals'])
        return Struct(point_reconstructed_ts=point_reconstructed_ts,
                      nf=nf, ts_modeintg_n=ts_modeintg_n,
                      crosssec_field=test_crossec_field)

analytical_dt = WGT.dt                  # Should be 0.0625

#orders = (1,2,3,4)
orders = (3,)
#dt_divisors = (128,64,32,16,8,4,2,1)
dt_divisors = (32,)
resses = {}



for order in orders:
    if not order in resses: resses[order]={}
    E_order = (order, True)
    B_order = (order-1, False) if order > 1 else (order, True)
    TRS = TestRun(mesh, BCs=Struct(E=freeE, B=freeB),
                  disc_orders=dict(E=E_order, B=B_order))
    TRS.setupSource()
    for dt_div in dt_divisors:
        print 'order: ', order
        dt = analytical_dt/dt_div
        no_steps = int(N.ceil(runtime/dt))
        #TRS.log_divisor = dt_div
        TRS.log_divisor = 1.
        TRS.set_dt(dt)      
        TRS.runSteps(no_steps)
        resses[order][dt_div] = TRS.getResult()


pickle.dump(resses, file('+coupledb_wg_direch_3rd_crossection.pickle', 'w'))

