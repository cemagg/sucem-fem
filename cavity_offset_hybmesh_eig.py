"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
from itertools import izip
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
from NewCode import Utilities, DifferentialForm
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Meshes.MeshIO import Femmesh
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import SubDimMesh, DifferentialForm, Waveforms, PostProc, Feeds, Runners
from NewCode.DifferentialForm import Discretiser, SubDimDiscretiser, \
     SubDimDiscretiserEntities, allconstrained
from NewCode.DiscretisedSystem import CurlCurlNewmarkDirechlet, BrickCurlCurlNewmarkDirechlet 
from NewCode.GeomGen.Hybrid import OffsetRect

from NewCode.Feeds import WaveguideEigenMatcher

import wg_hybrid_disc_fake
g_eps = 1e-10                           # Geometrical tollerance
eMAGUSImport.init('workspace')
tet_mesh = Mesh.Mesh(eMAGUSImport.get_listmesh())

print 'Tet_mesh elements: ', len(tet_mesh.elements)

# h0 = 1.0001/1.
# ha = 1.0001/4.
a,b,c, = geom_size = N.array([29,23,19], N.float64)
h0 = 1/1.
h = geom_size*h0
#tet_geom_size_q = N.array([1,1,1], N.int32)*2
#tet_geom_size_q = N.array([1,1,1], N.int32)
#tet_geom_size_q = N.array([1/h0/2,1/h0,1/h0], N.int32)
h[0] /= 2 ; tet_geom_size_q = N.array([1,1,1], N.int32)
volsize_q = N.array([2,1,1], N.int32)
#tet_offset_q = [1,1,1]
tet_offset_q = [0,0,0]
#volsize_q = N.int32(N.round(1/h0))*N.array([1,1,1], N.int32)

geom = OffsetRect(tet_geom_size_q, volsize_q, h, tet_offset_q)
geom.init_background_mesh()
geom.init_dead_background_element_set()
geom.init_geom()

brick_mesh = geom.background_mesh

order = 3
free = DifferentialForm.allfree

a_p, b_p, c_p = (close_to_point(xx, g_eps) for xx in (a,b,c))
zero_p = close_to_point(0, g_eps)

cb = DifferentialForm.constrained_on_boundary

def const_cav_bdry(ent):
    x,y,z = ent.nodeCoords.T
    return not (N.all(a_p(x)) or N.all(zero_p(x)) or 
                N.all(b_p(y)) or N.all(zero_p(y)) or
                N.all(c_p(z)) or N.all(zero_p(z)))

freeE = free

from analytical_WG_driver import WGT

drv_fun = WGT.discrete_drv_fn
z_measure = WGT.test_z
z_measure_p = close_to_point(z_measure, g_eps)
runtime = WGT.t_final
z_port = 0.
z_port_p = close_to_point(z_port, g_eps)
zero_p = close_to_point(0, g_eps)
on_port = lambda ent: N.all(z_port_p(ent.nodeCoords[:,2]))
def port_free(ent):
    x,y,z = ent.nodeCoords.T
    return not (N.all(a_p(x)) or N.all(zero_p(x)) or 
                N.all(b_p(y)) or N.all(zero_p(y)))
analytical_dt = WGT.dt

bdry_tri_geom = Femmesh.get_femmesh_tris(
    file('workspace/'+tet_mesh.FemmeshFilename))
bdry_ent_nodes = set([tuple(nds) for nds, mat in
                     izip(bdry_tri_geom.nodes, bdry_tri_geom.material)
                     if mat == 1])
bdry_ent_nodes.update(set([tuple(e_nds) for t_nds in bdry_ent_nodes
                           for e_nds in [ [t_nds[0], t_nds[1]],
                                          [t_nds[0], t_nds[2]],
                                          [t_nds[1], t_nds[2]] ]]))

hex_on_hbdry = geom.on_hbdry_hex
tet_on_hbdry = lambda ent: tuple(ent.nodes) in bdry_ent_nodes
on_hbdry = lambda ent: tet_on_hbdry(ent) if ent.meshtype == 'tet' \
           else hex_on_hbdry(ent)

hybrid_boundary_p = close_to_point(a/2, g_eps)
on_hbdry = lambda ent: N.all([hybrid_boundary_p(x) for (x,y,z) in ent.nodeCoords])

hyb_mesh = Struct(tet=tet_mesh, brick=brick_mesh, on_hbdry=on_hbdry)

class HybridCurlCurlNewmarkDirechlet(CurlCurlNewmarkDirechlet):
    DiscretiserModule = wg_hybrid_disc_fake
    DirechClass = BrickCurlCurlNewmarkDirechlet

    def setDirechBCs(self, DirechBCs, DirechVolBCs=allconstrained):
        self.direchSys = self.DirechClass(
            self.mesh.brick, order=self.order, BC=DirechBCs,
            volBC=DirechVolBCs, useQ=self.hasQ)
        self.direchSys.disc.setIntegrationRule(self.order*2-1)
        
class TestRun(Runners.TestRun):
    drv_fun = staticmethod(drv_fun)
    useLU = True
    SystemClass = HybridCurlCurlNewmarkDirechlet
    def __init__(self, mesh, **kwargs):
        self.system = self.SystemClass(mesh, **kwargs)

    def setupSource(self):
        self.system.drive_fun = self.drv_fun
        self.system.setDirechBCs(direch_free)
        self.sm = sm = Feeds.BrickSurfaceFieldMatcher()
        sm.initSubdim(self.system.disc.discs.brick, on_port, port_free)
        self.E_wg = E_wg = Feeds.gen_E_TE01(a,b)
        self.direch_dofArray = sm.matchKnown(E_wg)
        self.system.direchSys.dofs.dofArray[:] = self.direch_dofArray

    def setupLogging(self):
        divisor = self.log_divisor
        from NewCode import PostProc
        self.sm_m = sm_m = Feeds.BrickSurfaceFieldMatcher()
        sm_m.initSubdim(
            self.system.disc.discs.brick, on_measurement_port, measurement_port_free)
        self.measure_dofArray = sm_m.matchKnown(self.E_wg)
        self.tif = tif = Feeds.BrickTangentialFuncProjSurfaceIntegral(
            self.system.disc.discs.brick, on_measurement_port, measurement_port_free)
        self.test_pts = test_pts = N.array([[a/2, b/2, z_measure]], N.float64)
        test_elnos, test_el_coords = PostProc.LocatePoints(
            self.system.disc.discs.brick.mesh, test_pts)
        self.test_elnos = test_elnos ; self.test_el_coords = test_el_coords
        #self.system.addReconstructedLogger(test_elnos, test_el_coords, divisor=divisor)
        self.logged_dofnos = logged_dofnos = tif.superDOFMap
        self.system.addLogger(dofnos=logged_dofnos, divisor=divisor)
        
    def getResult(self):
        #rsv = N.array(self.system.loggedReconstructed[0]['vals'])[:,0,:]
        #point_reconstructed_ts = rsv[:,1]
        sm_m = self.sm_m
        measure_dofArray = self.measure_dofArray
        nf =  N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(measure_dofArray))
        ts_modeintg_n = N.array([
            N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(dof_vec))
            for dof_vec in self.system.loggedDOFs[tuple(self.logged_dofnos)].vals],
                                N.float64)/nf        
        return Struct(nf=nf, ts_modeintg_n=ts_modeintg_n)

analytical_dt = WGT.dt                  # Should be 0.0625

TRS = TestRun(hyb_mesh, order=order, BC=freeE, useQ=False,
              dead_elementset=geom.dead_background_elements)
M = TRS.system.disc.matrix.mass()
S = TRS.system.disc.matrix.stiffness()
st = TRS.system.disc.stuff
#M = st.M_brick_l
#S = st.S_brick_l

# hc_T_tetonly = st.hc_T[[i for i in range(st.hc_T.shape[0])
#                         if st.hc_T[i,:].nnz > 0]].copy()

# M = hc_T_tetonly.matmat(st.M_mtet.matmat(hc_T_tetonly.T))
# S = hc_T_tetonly.matmat(st.S_mtet.matmat(hc_T_tetonly.T))


from scipy.sparse.linalg.eigen.arpack import speigs
sigma = 0.01
def sigma_solve_gen(S):
    from NewCode.MatrixUtils import MatrixSolver
    ms = MatrixSolver(5e-9)
    cnt = [0]
    def solve(b):
        print 'Starting iterative solve'
        x = ms.solve_mat_vec(S, b)
        cnt[0] += 1
        print 'Done with iterative solve %i' %cnt[0]
        return x
    return solve
#sigma_solve = sigma_solve_gen(S - sigma*M)
print "(nodofs, nnz, sparsity %)", M.shape[0], M.nnz, M.nnz/M.shape[0]**2*100
#print 'Sparse LU decomposition'
#sigma_solve = scipy.sparse.linalg.dsolve.factorized(S - sigma*M)

print 'Solving Eigenproblem'
#w,v = speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, M.shape[0], sigma, 51, ncv=91)
w,v = scipy.linalg.eig(S.todense(), M.todense())
from AnalyticResults import PEC_cavity, accoustic_eigs, err_percentage

ares = PEC_cavity['rect-white']
#ares = accoustic_eigs['rect1x0.75x0.5']
#ares = PEC_cavity['rect1x0.25x5.1']

res =  N.array(sorted(N.abs(w[w > 0.0000001]))[0:10])

print res

err = err_percentage(ares, sorted(w[w > 0.0000001])[0:10])
RMS_err = Utilities.RMS(err[0:4])
#edges = eigsys.disc.geomEntities['edge']
#faces = eigsys.disc.geomEntities['face']

print RMS_err
print err
print 'Min eigvalue: ', N.min(w)
