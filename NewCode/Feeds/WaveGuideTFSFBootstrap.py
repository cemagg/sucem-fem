from __future__ import division

import numpy as N
import scipy
#
# Local Imports
#
import NewCode
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.Utilities import Struct, close_to_point
from NewCode import DifferentialForm, Waveforms,  Feeds
from NewCode.Feeds import BrickSurfaceFieldMatcher,\
     BrickTangentialFuncProjSurfaceIntegral
from NewCode.DiscretisedSystem.CoupledBonding import \
     CoupledBondingSystem, PMLBondedSystem, CoupledBondedSystem
from NewCode import PML

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
constrained = DifferentialForm.allconstrained

class WaveGuideTFSFBootstrap(object):
    g_eps = 1e-10                           # Geometrical tollerance
    useLU = True
    BondingSystemClass = CoupledBondingSystem
    DrivingSystemClass = CoupledBondedSystem
    PMLBondedSystemClass = PMLBondedSystem

    def __init__(self, a, b, h, order, drv_fun, no_PML_cells=20):
        self.a = a ; self.b = b; self.h = h
        self.order = order; self.drv_fun = drv_fun
        self.no_PML_cells = no_PML_cells
        self._calc_parms()
        self._calc_geomfuncs()
        self._setup_freefuns()
        E_order = (order, True)
        B_order = (order, True)
        self.disc_orders=dict(E=E_order, B=B_order)
        
    def _calc_parms(self):
        a,b,h = self.a, self.b, self.h
        no_PML_cells = self.no_PML_cells
        self.z_measure = z_measure = 3*h

        self.z_PML_goal = z_measure + h
        self.z_PML = N.ceil(self.z_PML_goal/h)*h
        self.c = (N.ceil(self.z_PML/h) + no_PML_cells)*h
        self.PML_len = no_PML_cells*h
        self.PML_sigma_max = PML.calc_PML_sigma_max(h)
        self.z_port = 0.
        self.drv_bond_surf_z = h
        
    def _calc_geomfuncs(self):
        g_eps = self.g_eps
        a,b,h = self.a, self.b, self.h
        no_PML_cells = self.no_PML_cells
        z_measure = self.z_measure
        self.z_port_p = close_to_point(self.z_port, g_eps)        
        self.z_measure_p = close_to_point(z_measure, g_eps)
        self.left_of_measure = lambda r: r[2] < z_measure
        self.z_PML_p = close_to_point(self.z_PML, g_eps)
        self.sigma_z_fn = PML.fn_1D2fn_3D(
            PML.PML_sigma_fn_1DR(self.PML_sigma_max,
                                 self.z_PML, self.PML_len), 2)
        self.sigma_fns = Struct(x=lambda r:0, y=lambda r:0, z=self.sigma_z_fn)
        self.zero_p = close_to_point(0, g_eps)
        self.a_p, self.b_p = (close_to_point(xx, g_eps) for xx in (a,b))
        self.on_drv_bond_surf_p = close_to_point(self.drv_bond_surf_z, g_eps)

    def in_bonding_sys(self, ent):
        return not (self.in_PML_sys(ent) or self.in_drv_sys(ent))

    def on_port(self, ent):
        return N.all(self.z_port_p(ent.nodeCoords[:,2]))

    def on_drv_bond_surf(self, ent):
        return N.all(self.on_drv_bond_surf_p(ent.nodeCoords[:,2]))

    def in_drv_sys(self, ent):
        return N.all(ent.nodeCoords[:,2] <= self.drv_bond_surf_z)

    def in_PML_sys(self, ent):
        return N.all(ent.nodeCoords[:,2] >= self.z_PML)

    def on_z_PML(self, ent):
        return N.all(self.z_PML_p(ent.nodeCoords[:,2]))

    def on_PML_bond_surf(self, ent):
        return self.on_z_PML(ent)
        
    def port_free(self, ent):
        ncx = ent.nodeCoords[:,0] ; ncy = ent.nodeCoords[:,1]
        return not (N.all(self.a_p(ncx)) or N.all(self.zero_p(ncx))
                    or N.all(self.b_p(ncy)) or N.all(self.zero_p(ncy)))

    def _setup_freefuns(self):
        self.freeE = freeE = cb
        self.freeB = lambda ent: freeE(ent) or self.on_port(ent)

    def on_measurement_port(self, ent):
        return N.all(self.z_measure_p(ent.nodeCoords[:,2]))

    def direch_free(self, ent):
        return self.on_port(ent) and self.port_free(ent)  

    def drv_sys_freeE(self, ent):
        return self.freeE(ent) and self.in_drv_sys(ent)

    def drv_sys_freeB(self, ent):
        return  self.freeB(ent) and self.in_drv_sys(ent)

    def pml_sys_free(self, ent):
        return self.freeE(ent) and self.in_PML_sys(ent)

    def bonding_sys_free(self, ent):
        return self.in_bonding_sys(ent) and self.freeE(ent)

    def init_mesh(self):
        self.mesh = BrickMesh.Mesh(
            BrickMeshGen.make_rect_cavity_brick_listmesh(
            self.a, self.b, self.c, self.h))
        print 'Mesh elements: ', len(self.mesh.elements)

    def init_systems(self):
        order = self.order
        print 'oooooooooorder: ', order
        in_bonding_sys = self.in_bonding_sys
        in_drv_sys = self.in_drv_sys
        mesh = self.mesh
        
        self.bonding_system = self.BondingSystemClass(in_bonding_sys)
        self.drv_system = self.DrivingSystemClass(
            mesh, BCs=Struct(E=self.drv_sys_freeE, B=self.drv_sys_freeB),
            volBCs=Struct(E=in_drv_sys, B=in_drv_sys),
            disc_orders=self.disc_orders, el_in_disc=in_drv_sys)
        self.bonding_system.add_bonded_system(self.drv_system, self.on_drv_bond_surf)
        self.pml_system = self.PMLBondedSystemClass(
            mesh, BCs=Struct(E=self.pml_sys_free, B=self.pml_sys_free),
            volBCs=Struct(E=self.in_PML_sys, B=self.in_PML_sys),
            disc_orders=self.disc_orders, el_in_disc=self.in_PML_sys)
        self.bonding_system.add_bonded_system(self.pml_system, self.on_PML_bond_surf)
        self.bonding_system.init_discs(
            mesh, BCs=Struct(E=self.bonding_sys_free, B=self.bonding_sys_free),
            volBCs=Struct(E=self.in_bonding_sys, B=self.in_bonding_sys),
            disc_orders=self.disc_orders, el_in_disc=self.in_bonding_sys)
        self.pml_system.set_bond(self.bonding_system.get_bond(self.pml_system))
        self.drv_system.set_bond(self.bonding_system.get_bond(self.drv_system))
        self.pml_system.set_sigmas(self.sigma_fns)
        self.systems = (self.bonding_system, self.drv_system, self.pml_system)

    def setupSource(self):
        order = self.order
        self.drv_system.drive_fun = self.drv_fun
        self.sm = sm = Feeds.BrickSurfaceFieldMatcher()
        sm.initSubdim(self.drv_system.discs.E, self.on_port, self.port_free)
        self.E_wg = E_wg = Feeds.gen_E_TE01(self.a, self.b)
        self.direch_dofArray = sm.matchKnown(E_wg)
        self.drv_system.setDirechBCs(Struct(E=self.direch_free, B=constrained))
        self.drv_system.direchSys.dofs.E.dofArray[:] = self.direch_dofArray
        
    def setupLogging(self):
        print "setupLogging!"
        self.logsys = logsys = self.bonding_system
        self.inc_field = Feeds.TFSFIncidentField(
            logsys.discs, self.on_measurement_port, self.left_of_measure)
        self.inc_field.init_incsurf(self.port_free)
        self.inc_field_info = self.inc_field.get_inc_info()
        self.log_dofnos = Struct(E=self.inc_field_info.E.dofnos,
                                 B=self.inc_field_info.B.dofnos)
        self.sm_m = sm_m = Feeds.BrickSurfaceFieldMatcher()
        sm_m.initSubdim(
            self.logsys.discs.E, self.on_measurement_port, self.port_free)
        self.measure_dofArray = sm_m.matchKnown(self.E_wg)
        self.nf = N.dot(self.measure_dofArray,
                        sm_m.subDisc.matrix.mass().matvec(self.measure_dofArray))
        self.log_dofarrays = Struct(E=self.logsys.dofs.E.dofArray,
                                    B=self.logsys.dofs.B.dofArray)

    def set_dt(self, dt):
        self.dt = dt
        self.resetHistory()
        for system in self.systems:
            if self.useLU: system.useLU = True
            print "Setting timestep"
            system.setTimestep(dt)
        print "Done setting timesteps"


    def init_steppers(self, runtime, yield_zero=True):
        self.runtime = runtime
        no_steps = int(N.ceil(runtime/self.dt))
        self.B_steppers = [system.step_B(no_steps) for system in self.systems]
        self.E_steppers = [system.step_E(no_steps) for system in self.systems]
        self.yield_zero = yield_zero # Return zero vals after runtime if true,
                                        # else raise StopIteration
    def next_drive(self):
        try:
            for bstep in self.B_steppers: bstep.next()
            for estep in self.E_steppers: estep.next()
            return self.extract_drive_values()
        except StopIteration:
            if self.yield_zero: return self.extract_drive_values(zero=True)
            else: raise
        

    def extract_drive_values(self, zero=False):
        E_dofs = self.log_dofarrays.E[self.log_dofnos.E]
        B_dofs = self.log_dofarrays.B[self.log_dofnos.B]
        if not zero:
            ts_modeintg_n = N.dot(
                self.measure_dofArray,
                self.sm_m.subDisc.matrix.mass().matvec(E_dofs)
                )/self.nf
        if zero:
            E_dofs = N.zeros_like(E_dofs)
            B_dofs = N.zeros_like(B_dofs)
            ts_modeintg_n = 0.
        return Struct(E_incdofs=E_dofs, B_incdofs=B_dofs,
                      mode_weighted=ts_modeintg_n)
            
    def get_stepsCompleted(self):
        return self.drv_system.n

    def resetHistory(self):
        for system in self.systems:
            system.reset_history()
        self.setupLogging()
        for system in self.systems:
            system.log()                # To add an entry for n=0    


