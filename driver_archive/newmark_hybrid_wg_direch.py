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
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import DifferentialForm, Waveforms, PostProc, Runners, Feeds, SystemMatrix
from NewCode.Feeds import BrickSurfaceFieldMatcher,\
     BrickTangentialFuncProjSurfaceIntegral

from NewCode.ImplicitExplicit import NewmarkHybridSystem

h = 1/8.
a,b,c = 1, 0.25, N.ceil(5.1/h)*h

mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [h, h, h]))

print 'Mesh elements: ', len(mesh.elements)

    
cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

freeE = cb

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
direch_free = lambda ent: on_port(ent) and port_free(ent)  


# cd = N.array([1,0,0], N.float64)
# test_crossection = PostProc.MakeLine(0.01*cd, cd*a*0.99, 51) + [0,0,0.3]
# newmarkSystem.addReconstructedLogger(
#     *PostProc.LocatePoints(newmarkSystem.disc.mesh, test_crossection))


class TestRun(Runners.TestRun):
    drv_fun = staticmethod(drv_fun)
    useLU = True
    SystemClass = NewmarkHybridSystem.NewmarkHybridSystem
    def __init__(self, mesh, **kwargs):
        BC = kwargs['BC']
        del(kwargs['BC'])
        order = kwargs['order']
        del(kwargs['order'])
        self.system = self.SystemClass(mesh, **kwargs)
        z_exp = 0.5
        implicit_elements = set(i for i, el in enumerate(mesh.elements)
                                if N.average(el.nodeCoords[:,2]) < z_exp )
        self.system.init_elgroups(implicit_elements)
        self.system.init_group_freefuns(BC)
        self.system.init_discs(order=order)
        self.system.init_block_matrices()
        self.system.init_merged_mats()
        self.system.init_dofs()        

    def setupSource(self):
        self.system.set_driveDOFs([0], 0., self.drv_fun)
        self.system.set_direchBCs(direch_free)
        self.sm = sm = BrickSurfaceFieldMatcher()
        sm.initSubdim(self.system.discs.E.a, on_port, port_free)
        self.E_wg = E_wg = Feeds.gen_E_TE01(a,b)
        self.direch_dofArray = sm.matchKnown(E_wg)
        self.system.direch_dofs.dofArray[:] = self.direch_dofArray

    def setupLogging(self):
        divisor = self.log_divisor
        from NewCode import PostProc
        self.sm_m = sm_m = BrickSurfaceFieldMatcher()
        sm_m.initSubdim(self.system.discs.E.e, on_measurement_port, measurement_port_free)
        self.measure_dofArray = sm_m.matchKnown(self.E_wg)
        self.tif = tif = BrickTangentialFuncProjSurfaceIntegral(
            self.system.discs.E.e, on_measurement_port, measurement_port_free)
        self.test_pts = test_pts = N.array([[a/2, b/2, z_measure]], N.float64)
        test_elnos, test_el_coords = PostProc.LocatePoints(self.system.discs.E.a.mesh, test_pts)
        self.test_elnos = test_elnos ; self.test_el_coords = test_el_coords
        self.system.addReconstructedLogger('e', test_elnos, test_el_coords, divisor=divisor)
        self.logged_dofnos = logged_dofnos = tif.superDOFMap
        self.system.addLogger(group='e', dofnos=logged_dofnos, divisor=divisor)
        
    def getResult(self):
        rsv = N.array(self.system.loggedReconstructed[0]['vals'])[:,0,:]
        point_reconstructed_ts = rsv[:,1]
        sm_m = self.sm_m
        measure_dofArray = self.measure_dofArray
        nf =  N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(measure_dofArray))
        ts_modeintg_n = N.array([
            N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(dof_vec))
            for dof_vec in self.system.loggedDOFs['e', tuple(self.logged_dofnos)].vals],
                                N.float64)/nf
        return Struct(point_reconstructed_ts=point_reconstructed_ts,
                      nf=nf, ts_modeintg_n=ts_modeintg_n)

analytical_dt = WGT.dt                  # Should be 0.0625

#orders = (4,3,2,1)
orders = (2,1,)
#dt_divisors = (128,64,32,16,8,4,2,1)
dt_divisors = {1:(128,2,3,1),
               2:(128,3),
               3:(128,4),
               4:(128,64,32,16,8,6)}
resses = {}
for order in orders:
    if not order in resses: resses[order]={}
    TRS = TestRun(mesh, BC=freeE, implicit_beta=0.25001, order=order)
    TRS.setupSource()    
    for dt_div in dt_divisors[order]:
        dt = analytical_dt/dt_div
        no_steps = int(N.ceil(runtime/dt))
        TRS.log_divisor = dt_div
        TRS.set_dt(dt)      
        TRS.runSteps(no_steps)
        resses[order][dt_div] = TRS.getResult()


pickle.dump(resses, file('+newmark_hybrid_tmp.pickle', 'w'))

# el0 = TRS.system.discs.E.a.elements[0]
# el0_D = TRS.system.discs.E.a.D().elements[0]

# M_e = SystemMatrix.local_self_projection_matrix(el0)
# S_e = SystemMatrix.local_self_projection_matrix(el0_D)
# w_max = N.max(N.abs(scipy.linalg.eigvals(S_e,M_e)))
# dt_max = 2/N.sqrt(w_max)
