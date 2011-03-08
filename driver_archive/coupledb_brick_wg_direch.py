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
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import DifferentialForm, Waveforms, PostProc, Runners, Feeds, SystemMatrix
from NewCode.Feeds import BrickSurfaceFieldMatcher,\
     BrickTangentialFuncProjSurfaceIntegral
from NewCode.DiscretisedSystem import BrickCoupledFirstOrderSystemBDirechlet

h = 1/8.
wg_len = 5.1
a,b,c = 1, 0.25, N.ceil(wg_len/h)*h

mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [h, h, h]))

print 'Mesh elements: ', len(mesh.elements)

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
def port_free(ent):
    ncx = ent.nodeCoords[:,0] ; ncy = ent.nodeCoords[:,1]
    return not (N.all(a_p(ncx)) or N.all(zero_p(ncx))
                or N.all(b_p(ncy)) or N.all(zero_p(ncy)))
                                        
on_measurement_port = lambda ent: N.all(z_measure_p(ent.nodeCoords[:,2]))
measurement_port_free = port_free
freeE = cb
freeB = lambda ent: freeE(ent) or on_port(ent)
direch_free = lambda ent: on_port(ent) and port_free(ent)  

class TestRun(Runners.TestRun):
    drv_fun = staticmethod(drv_fun)
    useLU = True
    SystemClass = BrickCoupledFirstOrderSystemBDirechlet

    def __init__(self, mesh, **kwargs):
        order = kwargs['disc_orders']['E'][0]
        print 'oooooooooorder: ', order
        self.system = self.SystemClass(mesh, **kwargs)
        #self.system.discs.E.diagonalise()
        #self.system.discs.B.diagonalise()
        
    def setupSource(self):
        self.system.drive_fun = self.drv_fun
        self.sm = sm = Feeds.BrickSurfaceFieldMatcher()
        sm.initSubdim(self.system.discs.E, on_port, port_free)
        self.E_wg = E_wg = Feeds.gen_E_TE01(a,b)
        self.direch_dofArray = sm.matchKnown(E_wg)
        self.system.setDirechBCs(Struct(E=direch_free, B=constrained))
        #self.system.direchSys.discs.E.diagonalise()
        self.system.direchSys.dofs.E.dofArray[:] = self.direch_dofArray

    def setupLogging(self):
        divisor = self.log_divisor
        from NewCode import PostProc
        self.sm_m = sm_m = Feeds.BrickSurfaceFieldMatcher()
        sm_m.initSubdim(self.system.discs.E, on_measurement_port, measurement_port_free)
        self.measure_dofArray = sm_m.matchKnown(self.E_wg)
        self.tif = tif = Feeds.BrickTangentialFuncProjSurfaceIntegral(
            self.system.discs.E, on_measurement_port, measurement_port_free)
        self.test_pts = test_pts = N.array([[a/2, b/2, z_measure]], N.float64)
        test_elnos, test_el_coords = PostProc.LocatePoints(self.system.discs.E.mesh, test_pts)
        self.test_elnos = test_elnos ; self.test_el_coords = test_el_coords
        self.system.addReconstructedLogger('E', test_elnos, test_el_coords, divisor=divisor)
        self.logged_dofnos = logged_dofnos = tif.superDOFMap
        self.system.addLogger('E', dofnos=logged_dofnos, divisor=divisor)

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
        return Struct(point_reconstructed_ts=point_reconstructed_ts,
                      nf=nf, ts_modeintg_n=ts_modeintg_n)

analytical_dt = WGT.dt                  # Should be 0.0625

#orders = (1,2,3,4)
orders = (3,2,1)
#orders = (1,)
# For consistent mass
dt_divisors = {1:(128,2),
               2:(128,4),
               3:(128,6),
               4:(128,9)}
# For lumped mass h = 1/8 
# dt_divisors = {1:(128,1),
#                2:(128,3),
#                3:(128,4),
#                4:(128,6)}
resses = {}
for order in orders:
    if not order in resses: resses[order]={}
    E_order = (order, True)
    B_order = (order, True)
    TRS = TestRun(mesh, BCs=Struct(E=freeE, B=freeB),
                  disc_orders=dict(E=E_order, B=B_order),
                  btype='cohen98')
    TRS.setupSource()
    for dt_div in dt_divisors[order]:
        print 'order: ', order
        dt = analytical_dt/dt_div
        no_steps = int(N.ceil(runtime/dt))
        TRS.log_divisor = dt_div
        TRS.set_dt(dt)      
        TRS.runSteps(no_steps)
        resses[order][dt_div] = TRS.getResult()


pickle.dump(resses, file('+coupledb_brick_wg_tmp.pickle', 'w'))

