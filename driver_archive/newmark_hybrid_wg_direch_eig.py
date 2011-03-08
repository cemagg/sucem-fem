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
from NewCode import Utilities
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import DifferentialForm, Waveforms, PostProc, Runners, Feeds, SystemMatrix
from NewCode.Feeds import BrickSurfaceFieldMatcher,\
     BrickTangentialFuncProjSurfaceIntegral

from NewCode.ImplicitExplicit import NewmarkHybridSystem

h = 1/2.
a,b,c = 1, 0.25, 5.1

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
        z_exp = N.ceil(0.5/mesh.gridStepSize[2])*mesh.gridStepSize[2]
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

order = 3
TRS = TestRun(mesh, BC=freeE, implicit_beta=0, order=order)
TRS.system.set_dt(1)      
M = TRS.system.merged_matrices.A()
S = 2*M - TRS.system.merged_matrices.B()


from scipy.sparse.linalg.eigen.arpack import speigs
sigma = 0.1
print "(nodofs, nnz, sparsity %)", M.shape[0], M.nnz, M.nnz/M.shape[0]**2.0*100
print 'Sparse LU decomposition'
sigma_solve = scipy.sparse.linalg.dsolve.factorized(S - sigma*M)

print 'Solving Eigenproblem'
#w,v = speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, M.shape[0], sigma, 51, ncv=91)
w,v = scipy.linalg.eig(S.todense(), M.todense())

res =  N.array(sorted(N.abs(w[w > 0.0000001]))[0:10])

from AnalyticResults import PEC_cavity, accoustic_eigs, err_percentage
ares = PEC_cavity['rect1x0.25x5.1']

err = err_percentage(ares, res)
RMS_err = Utilities.RMS(err)

# err_mtet = err_percentage(ares, res_mtet)
# RMS_err_mtet = Utilities.RMS(err_mtet)

print RMS_err
print err
print 'Min real eigvalue: ', N.min(N.real(w))
