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
from NewCode.DiscretisedSystem.CoupledBonding import \
     CoupledBondingSystem, PMLBondedSystem, CoupledBondedSystem
from analytical_WG_driver import WGT
g_eps = 1e-10                           # Geometrical tollerance
drv_fun = WGT.discrete_drv_fn
z_measure = WGT.test_z
z_measure_p = close_to_point(z_measure, g_eps)

h = 1/4.
z_PML_goal = z_measure
no_PML_cells = 10
PML_m = 3
sigma_factor = 1
#a,b,c = 1, 0.25, N.ceil(5.1/h)*h
a,b,c = 1., 0.25, (N.ceil(z_PML_goal/h) + no_PML_cells)*h
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
constrained = DifferentialForm.allconstrained
 

runtime = WGT.t_final
z_port = 0.
z_port_p = close_to_point(z_port, g_eps)
zero_p = close_to_point(0, g_eps)
a_p, b_p = (close_to_point(xx, g_eps) for xx in (a,b))
on_port = lambda ent: N.all(z_port_p(ent.nodeCoords[:,2]))

in_drv_sys = lambda ent: N.all(ent.nodeCoords[:,2] <= h)
on_drv_bond_surf_p = close_to_point(h, g_eps)
on_drv_bond_surf = lambda ent: N.all(on_drv_bond_surf_p(ent.nodeCoords[:,2]))
in_PML_sys = lambda ent: N.all(ent.nodeCoords[:,2] >= z_PML)
on_z_PML = lambda ent: N.all(z_PML_p(ent.nodeCoords[:,2]))
on_PML_bond_surf = on_z_PML                               
in_bonding_sys = lambda ent: not (in_PML_sys(ent) or in_drv_sys(ent))

def port_free(ent):
    ncx = ent.nodeCoords[:,0] ; ncy = ent.nodeCoords[:,1]
    return not (N.all(a_p(ncx)) or N.all(zero_p(ncx))
                or N.all(b_p(ncy)) or N.all(zero_p(ncy)))
                                        
on_measurement_port = lambda ent: N.all(z_measure_p(ent.nodeCoords[:,2]))
measurement_port_free = port_free
freeE = cb
freeB = lambda ent: freeE(ent) or on_port(ent)
direch_free = lambda ent: on_port(ent) and port_free(ent)  

drv_sys_freeE = lambda ent: freeE(ent) and in_drv_sys(ent)
drv_sys_freeB = lambda ent: freeB(ent) and in_drv_sys(ent)
pml_sys_free = lambda ent: freeE(ent) and in_PML_sys(ent)
bonding_sys_free = lambda ent: in_bonding_sys(ent) and freeE(ent)

class TestRun(Runners.TestRun):
    drv_fun = staticmethod(drv_fun)
    useLU = True
    BondingSystemClass = CoupledBondingSystem
    DrivingSystemClass = CoupledBondedSystem
    PMLBondedSystemClass = PMLBondedSystem
    sigma_fns = Struct(x=lambda r:0, y=lambda r:0, z=sigma_z_fn)
    
    def __init__(self, mesh, **kwargs):
        self.order = order = kwargs['disc_orders']['E'][0]
        print 'oooooooooorder: ', order
        self.bonding_system = self.BondingSystemClass(in_bonding_sys)
        self.drv_system = self.DrivingSystemClass(
            mesh, BCs=Struct(E=drv_sys_freeE, B=drv_sys_freeB),
            volBCs=Struct(E=in_drv_sys, B=in_drv_sys),
            disc_orders=kwargs['disc_orders'], el_in_disc=in_drv_sys)
        self.bonding_system.add_bonded_system(self.drv_system, on_drv_bond_surf)
        self.pml_system = self.PMLBondedSystemClass(
            mesh, BCs=Struct(E=pml_sys_free, B=pml_sys_free),
            volBCs=Struct(E=in_PML_sys, B=in_PML_sys),
            disc_orders=kwargs['disc_orders'], el_in_disc=in_PML_sys)
        self.bonding_system.add_bonded_system(self.pml_system, on_PML_bond_surf)
        self.bonding_system.init_discs(
            mesh, BCs=Struct(E=bonding_sys_free, B=bonding_sys_free),
            volBCs=Struct(E=in_bonding_sys, B=in_bonding_sys),
            disc_orders=kwargs['disc_orders'], el_in_disc=in_bonding_sys)
        self.pml_system.set_bond(self.bonding_system.get_bond(self.pml_system))
        self.drv_system.set_bond(self.bonding_system.get_bond(self.drv_system))
        self.setupPML()
        self.systems = (self.bonding_system, self.drv_system, self.pml_system)

    def setupPML(self):
        self.pml_system.set_sigmas(self.sigma_fns)

    def setupSource(self):
        order = self.order
        self.drv_system.drive_fun = self.drv_fun
        self.sm = sm = Feeds.BrickSurfaceFieldMatcher()
        sm.initSubdim(self.drv_system.discs.E, on_port, port_free)
        self.E_wg = E_wg = Feeds.gen_E_TE01(a,b)
        self.direch_dofArray = sm.matchKnown(E_wg)
        self.drv_system.setDirechBCs(Struct(E=direch_free, B=constrained))
        self.drv_system.direchSys.dofs.E.dofArray[:] = self.direch_dofArray

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
        return self.drv_system.n

    def resetHistory(self):
        for system in self.systems:
            system.reset_history()
        self.setupLogging()
        for system in self.systems:
            system.log()                # To add an entry for n=0    

analytical_dt = WGT.dt                  # Should be 0.0625

#orders = (1,2,3,4)
#orders = (3,2,1)
orders = (1,2)
# For lumped mass h = 1/8 (not sure about 4th order)
dt_divisors = {1:(1,),
               2:(3,),
               3:(4,),
               4:(128,9)}
resses = {}
for order in orders:
    if not order in resses: resses[order]={}
    E_order = (order, True)
    B_order = (order, True)
    TRS = TestRun(mesh, disc_orders=dict(E=E_order, B=B_order),
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


#pickle.dump(resses, file('+brick_4_wg_PML_tmp.pickle', 'w'))

