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
from numpy.testing import assert_equal, assert_almost_equal, assert_approx_equal

#
# Local Imports
#
import NewCode
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import DifferentialForm, Waveforms, PostProc, Runners, Feeds, SystemMatrix
from NewCode.Feeds import BrickSurfaceFieldMatcher,\
     BrickTangentialFuncProjSurfaceIntegral

from NewCode.ImplicitExplicit.HybridMeshNewmarkHybridSystem import \
     HybridMeshNewmarkHybridSystem

eMAGUSImport.init('workspace')
tet_mesh = Mesh.Mesh(eMAGUSImport.get_listmesh())

print 'Tet_mesh elements: ', len(tet_mesh.elements)

bfrac = 0.45
h = 1/4.
a,b,c = 1, 0.25, 5.1
bfrac = N.ceil(bfrac*c/h)*h/c

brick_mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c*bfrac, [h, h, h]))

print 'Brick-Mesh elements: ', len(brick_mesh.elements)
    
cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

g_eps = 1e-10                           # Geometrical tollerance
hybrid_boundary_p = close_to_point(c*bfrac, g_eps)
on_hbdry = lambda ent: N.all([hybrid_boundary_p(z) for (x,y,z) in ent.nodeCoords])
a_p, b_p, c_p = (close_to_point(xx, g_eps) for xx in (a,b,c))
zero_p = close_to_point(0, g_eps)

def freeE(ent):
    x,y,z = ent.nodeCoords.T
    return not (N.all(a_p(x)) or N.all(zero_p(x)) or 
                N.all(b_p(y)) or N.all(zero_p(y)) or
                N.all(c_p(z)) or N.all(zero_p(z)))

z_port = 0.
z_port_p = close_to_point(z_port, g_eps)
on_port = lambda ent: N.all(z_port_p(ent.nodeCoords[:,2]))

from analytical_WG_driver import WGT


drv_fun = WGT.discrete_drv_fn
z_measure = WGT.test_z
z_measure_p = close_to_point(z_measure, g_eps)
runtime = WGT.t_final
z_port = 0.
z_port_p = close_to_point(z_port, g_eps)
a = 1. ; b = 0.25                           # WG x/y dim

def port_free(ent):
    x,y,z = ent.nodeCoords.T
    return not (N.all(a_p(x)) or N.all(zero_p(x)) or 
                N.all(b_p(y)) or N.all(zero_p(y)))


on_measurement_port = lambda ent: N.all(z_measure_p(ent.nodeCoords[:,2]))
measurement_port_free = port_free
direch_free = lambda ent: on_port(ent) and port_free(ent)  


class TestRun(Runners.TestRun):
    drv_fun = staticmethod(drv_fun)
    useLU = True
    SystemClass = HybridMeshNewmarkHybridSystem
    drive_group = 'e'
    log_group = 'e'
    def __init__(self, tet_mesh, brick_mesh, BC, on_hbdry, order, **kwargs):
        self.system = self.SystemClass(tet_mesh, brick_mesh, **kwargs)
        self.system.init_group_freefuns(on_hbdry, BC)
        self.system.init_discs(order=order)
        self.system.init_block_matrices()
        self.system.init_merged_mats()
        self.system.init_dofs()        

    def setupSource(self):
        self.system.set_driveDOFs(self.drive_group, [0], 0., self.drv_fun)
        self.system.set_direchBCs(direch_free)
        self.sm = sm = BrickSurfaceFieldMatcher()
        sm.initSubdim(self.system.discs.E[self.drive_group], on_port, port_free)
        self.E_wg = E_wg = Feeds.gen_E_TE01(a,b)
        self.direch_dofArray = sm.matchKnown(E_wg)
        self.system.direch_dofs.dofArray[:] = self.direch_dofArray

    def setupLogging(self):
        divisor = self.log_divisor
        lg = self.log_group ; log_disc = self.system.discs.E[lg]
        from NewCode import PostProc
        self.sm_m = sm_m = BrickSurfaceFieldMatcher()
        sm_m.initSubdim(log_disc, on_measurement_port, measurement_port_free)
        self.measure_dofArray = sm_m.matchKnown(self.E_wg)
        self.tif = tif = BrickTangentialFuncProjSurfaceIntegral(
            log_disc, on_measurement_port, measurement_port_free)
        self.test_pts = test_pts = N.array([[a/2, b/2, z_measure]], N.float64)
        test_elnos, test_el_coords = PostProc.LocatePoints(log_disc.mesh, test_pts)
        self.test_elnos = test_elnos ; self.test_el_coords = test_el_coords
        self.system.addReconstructedLogger(lg, test_elnos, test_el_coords, divisor=divisor)
        self.logged_dofnos = logged_dofnos = tif.superDOFMap
        self.system.addLogger(group=lg, dofnos=logged_dofnos, divisor=divisor)
        
    def getResult(self):
        lg = self.log_group 
        rsv = N.array(self.system.loggedReconstructed[0]['vals'])[:,0,:]
        point_reconstructed_ts = rsv[:,1]
        sm_m = self.sm_m
        measure_dofArray = self.measure_dofArray
        nf =  N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(measure_dofArray))
        ts_modeintg_n = N.array([
            N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(dof_vec))
            for dof_vec in self.system.loggedDOFs[lg, tuple(self.logged_dofnos)].vals],
                                N.float64)/nf
        return Struct(point_reconstructed_ts=point_reconstructed_ts,
                      nf=nf, ts_modeintg_n=ts_modeintg_n)

analytical_dt = WGT.dt                  # Should be 0.0625

#orders = (4,3,2,1)
orders = (3,2,1)
#dt_divisors = (128,64,32,16,8,4,2,1)
dt_divisors = {1:(128,1),
               2:(128,3),
               3:(128,4),
               4:(128,64,32,16,8,6)}
#dt_divisors = (16,)
resses = {}
for order in orders:
    if not order in resses: resses[order]={}
    TRS = TestRun(tet_mesh, brick_mesh, BC=freeE, on_hbdry=on_hbdry,
                  implicit_beta=0.25, order=order)
    TRS.setupSource()    
    for dt_div in dt_divisors[order]:
        dt = analytical_dt/dt_div
        no_steps = int(N.ceil(runtime/dt))
        TRS.log_divisor = dt_div
        TRS.set_dt(dt)      
        TRS.runSteps(no_steps)
        resses[order][dt_div] = TRS.getResult()


pickle.dump(resses, file('+newmark_hybrid_wg_hybmesh-tmp.pickle', 'w'))

# el0 = TRS.system.discs.E.a.elements[0]
# el0_D = TRS.system.discs.E.a.D().elements[0]

# M_e = SystemMatrix.local_self_projection_matrix(el0)
# S_e = SystemMatrix.local_self_projection_matrix(el0_D)
# w_max = N.max(N.abs(scipy.linalg.eigvals(S_e,M_e)))
# dt_max = 2/N.sqrt(w_max)
