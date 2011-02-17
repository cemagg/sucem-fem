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
from NewCode import BrickSubDimMesh, DifferentialForm, Waveforms, PostProc, Feeds, Runners
from NewCode.ImplicitExplicit import NewmarkBondedSystem
from NewCode.DiscretisedSystem.CoupledBonding import \
     CoupledBondingSystem, PMLBondedSystem
from NewCode.DifferentialForm import BrickSubDimDiscretiser, BrickSubDimDiscretiserEntities
from NewCode.Feeds import BrickSurfaceFieldMatcher,\
     BrickTangentialFuncProjSurfaceIntegral

g_eps = 1e-10                           # Geometrical tollerance
from analytical_WG_driver import WGT
drv_fun = WGT.discrete_drv_fn
z_measure = WGT.test_z
z_measure_p = close_to_point(z_measure, g_eps)
runtime = WGT.t_final

h = 1/4.

z_PML_goal = z_measure
no_PML_cells = 10
PML_m = 3
sigma_factor = 1
#a,b,c = 1, 0.25, N.ceil(5.1/h)*h
a,b,c = 1, 0.25, (N.ceil(z_PML_goal/h) + no_PML_cells)*h
z_PML = N.ceil(z_PML_goal/h)*h ; PML_len = no_PML_cells*h
wave_imp = 1.
PML_sigma_max = 0.8*(PML_m+1)/h/wave_imp*sigma_factor
z_PML_p = close_to_point(z_PML, g_eps)

sigma_z_fn = lambda r: N.where(
    N.abs(r.T[2]) >= z_PML,
    PML_sigma_max*((N.abs(r.T[2])-z_PML)/PML_len)**PML_m,
    0.)

mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [h, h, h]))

print 'Mesh elements: ', len(mesh.elements)


cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

freeE = cb


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

in_PML_sys = lambda ent: N.all(ent.nodeCoords[:,2] >= z_PML)
on_z_PML = lambda ent: N.all(z_PML_p(ent.nodeCoords[:,2]))
on_PML_bond_surf = on_z_PML                               
pml_sys_free = lambda ent: in_PML_sys(ent) and freeE(ent)

class TestRun(Runners.TestRun):
    drv_fun = staticmethod(drv_fun)
    useLU = True
    HybridSystemClass = NewmarkBondedSystem.NewmarkBondedSystem
    BondingSystemClass = CoupledBondingSystem
    PMLBondedSystemClass = PMLBondedSystem
    sigma_fns = Struct(x=lambda r:0, y=lambda r:0, z=sigma_z_fn)
    def __init__(self, mesh, **kwargs):
        BC = kwargs['BC']
        del(kwargs['BC'])
        order = kwargs['order']
        orders = Struct(E=(order, True), B=(order, True))
        del(kwargs['order'])
        self.hybrid_system = self.HybridSystemClass(mesh, **kwargs)
        z_exp = 0.5
        implicit_elements = set(i for i, el in enumerate(mesh.elements)
                                if N.average(el.nodeCoords[:,2]) < z_exp )
        hybsys_bdry = z_exp + h
        on_hybsys_bdry_p = close_to_point(hybsys_bdry, g_eps)
        on_hybsys_bond_surf = lambda ent: N.all(
            on_hybsys_bdry_p(ent.nodeCoords[:,2]))
        in_bonding_sys = lambda ent: N.all(ent.nodeCoords[:,2] <= z_PML) \
                         and N.all(ent.nodeCoords[:,2] >= hybsys_bdry)
        self.hybrid_system.init_elgroups(implicit_elements)
        self.hybrid_system.init_group_freefuns(BC)
        self.hybrid_system.init_discs(order=order)
        self.hybrid_system.init_block_matrices()
        self.hybrid_system.init_merged_mats()
        self.hybrid_system.init_dofs()        
        self.pml_system = self.PMLBondedSystemClass(
            mesh, BCs=Struct(E=pml_sys_free, B=pml_sys_free),
            volBCs=Struct(E=in_PML_sys, B=in_PML_sys),
            disc_orders=orders, el_in_disc=in_PML_sys)        
        self.bonding_system = self.BondingSystemClass(in_bonding_sys)
        self.bonding_system.add_bonded_system(self.pml_system, on_PML_bond_surf)
        self.bonding_system.add_bonded_system(self.hybrid_system, on_hybsys_bond_surf)
        freeE_e = self.hybrid_system.group_freefuns.e[0]
        bonding_sys_free = lambda ent: in_bonding_sys(ent) and freeE_e(ent)
        self.bonding_system.init_discs(
            mesh, BCs=Struct(E=bonding_sys_free, B=bonding_sys_free),
            volBCs=Struct(E=in_bonding_sys, B=in_bonding_sys),
            disc_orders=orders, el_in_disc=in_bonding_sys)
        self.pml_system.set_bond(self.bonding_system.get_bond(self.pml_system))
        self.hybrid_system.set_bond(self.bonding_system.get_bond(self.hybrid_system))
        self.pml_system.set_sigmas(self.sigma_fns)
        self.systems = (self.hybrid_system, self.bonding_system, self.pml_system)
        
    def setupSource(self):
        self.hybrid_system.set_driveDOFs([0], 0., self.drv_fun)
        self.hybrid_system.set_direchBCs(direch_free)
        self.sm = sm = BrickSurfaceFieldMatcher()
        sm.initSubdim(self.hybrid_system.discs.E.a, on_port, port_free)
        self.E_wg = E_wg = Feeds.gen_E_TE01(a,b)
        self.direch_dofArray = sm.matchKnown(E_wg)
        self.hybrid_system.direch_dofs.dofArray[:] = self.direch_dofArray

    def setupLogging(self):
        divisor = self.log_divisor
        from NewCode import PostProc
        self.sm_m = sm_m = Feeds.BrickSurfaceFieldMatcher()
        sm_m.initSubdim(self.pml_system.discs.E, on_measurement_port, measurement_port_free)
        self.measure_dofArray = sm_m.matchKnown(self.E_wg)
        self.tif = tif = Feeds.BrickTangentialFuncProjSurfaceIntegral(
            self.pml_system.discs.E, on_measurement_port, measurement_port_free)
        self.logged_dofnos = logged_dofnos = tif.superDOFMap
        self.pml_system.addLogger('E', dofnos=logged_dofnos, divisor=divisor)

    def getResult(self):
        sm_m = self.sm_m
        measure_dofArray = self.measure_dofArray
        nf =  N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(measure_dofArray))
        ts_modeintg_n = N.array([
            N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(dof_vec))
            for dof_vec in self.pml_system.loggedDOFs.E[tuple(self.logged_dofnos)].vals],
                                N.float64)/nf        
        return Struct(nf=nf, ts_modeintg_n=ts_modeintg_n)

    def set_dt(self, dt):
        self.dt = dt
        self.resetHistory()
        for system in self.systems:
            if self.useLU: system.useLU = True
            print "Setting timestep"
            system.setTimestep(dt)
        print "Done setting timesteps"

    def runSteps(self, n_steps):
        B_steppers = [system.step_B(n_steps) for system in self.systems]
        E_steppers = [system.step_E(n_steps) for system in self.systems]
        for n in xrange(n_steps):
            for bstep in B_steppers: bstep.next()
            for estep in E_steppers: estep.next()
        
    def get_stepsCompleted(self):
        return self.hybrid_system.n

    def resetHistory(self):
        for system in self.systems:
            system.reset_history()
        self.setupLogging()
        for system in self.systems:
            system.log()                # To add an entry for n=0    
    

analytical_dt = WGT.dt                  # Should be 0.0625

orders = (3,2,1,)
# For lumped mass h = 1/8 
dt_divisors = {1:(128,1),
               2:(128,3),
               3:(128,4),
               4:(128,6)}
resses = {}
for order in orders:
    if not order in resses: resses[order]={}
    TRS = TestRun(mesh, BC=freeE, implicit_beta=0.25001, order=order)
    TRS.setupSource()    
    for dt_div in dt_divisors[order]:
        dt = analytical_dt/dt_div
        no_steps = int(N.ceil(runtime/dt))#*100
        TRS.log_divisor = dt_div
        TRS.set_dt(dt)      
        TRS.runSteps(no_steps)
        resses[order][dt_div] = TRS.getResult()

bsdE = TRS.bonding_system.discs.E
hsdEd = TRS.hybrid_system.discs.E.d
bsdB = TRS.bonding_system.discs.B
hsdBd = TRS.hybrid_system.discs.B.d
psdB = TRS.pml_system.discs.B
psdE = TRS.pml_system.discs.E

E_discs = (bsdE, psdE, hsdEd)
B_discs = (bsdB, psdB, hsdBd)

edge_free_tables = [disc.permuter.globalEntityPermutationTable('edge').T >= 0
                    for disc in E_discs]

face_free_tables = [disc.permuter.globalEntityPermutationTable('face').T >= 0
                    for disc in B_discs]

pickle.dump(resses, file('+newmark_pml_hybrid_tmp.pickle', 'w'))

