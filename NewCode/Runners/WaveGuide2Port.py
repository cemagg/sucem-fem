from __future__ import division
import numpy as N

from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.Utilities import Struct, partial, close_to_point, almost_leq, almost_geq
from NewCode.DiscretisedSystem.CoupledBonding import \
     CoupledBondingSystem, PMLBondedSystem, CoupledBondedSystem
from NewCode.Feeds.WaveGuideTFSFBootstrap import WaveGuideTFSFBootstrap
from NewCode.Feeds import BrickSurfaceFieldMatcher,\
     BrickTangentialFuncProjSurfaceIntegral
from NewCode import DifferentialForm,  Feeds
from NewCode import PML

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
constrained = DifferentialForm.allconstrained

class WaveGuide2Port(object):
    g_eps = 1e-10                       # Geometrical tollerance
    useLU = True
    BondingSystemClass = CoupledBondingSystem
    PMLBondedSystemClass = PMLBondedSystem
    source_no_PML_cells = 25
    source_runtime = 100000             # Essentially infinite
    short = False
    measure1_on_incport = False

    def __init__(self, a, b, h, order, drv_fun, no_PML_cells):
        self.a = a ; self.b = b; self.h = h
        self.order = order; self.drv_fun = drv_fun
        self.no_PML_cells = no_PML_cells
        E_order = (order, True)
        B_order = (order, True)
        self.disc_orders=dict(E=E_order, B=B_order)

    def set_dims(self, z_incport, z_measure1, z_measure2, short_after=None):
        """Note z_measure2 > z_incport >= z_measure1 should be true and the
        distance between z_measure1 and z_incport should be zero or a multiple
        of h"""
        a,b,h = self.a, self.b, self.h
        no_PML_cells = self.no_PML_cells
        assert(z_measure2 > z_incport)
        assert(z_measure1 <= z_incport)
        if N.abs(z_measure1 - z_incport) < self.g_eps:
            self.measure1_on_incport = True
        self.z_incport = z_incport
        self.z_measure1 = z_measure1
        self.z_measure2 = z_measure2
        self.z_PML_goal1 = self._calc_z_PML1_goal()
        self.z_PML1 = N.floor(self.z_PML_goal1/h)*h
        if short_after is None:
            self.z_PML_goal2 = self._calc_z_PML2_goal()
            self.z_PML2 = N.ceil(self.z_PML_goal2/h)*h
            self.c = (-self.z_PML1/h + self.z_PML2/h + no_PML_cells*2)*h
        else:
            self.thru_len = N.ceil(short_after/h)*h
            self.short = True
            self.z_short = self.z_incport + self.thru_len
            self.c = (-self.z_PML1/h + self.thru_len/h + no_PML_cells)*h
            self.z_PML2 = 1e99           # "Infinite"
        self.PML_len = no_PML_cells*h
        self.PML_sigma_max = PML.calc_PML_sigma_max(h)
        self.mesh_offset = self.z_PML1-self.PML_len
        self.z_guide_start = self.mesh_offset
        self.z_guide_end = self.c + self.mesh_offset
        self._calc_geomfuncs()
        self._setup_freefuns()

    def _calc_geomfuncs(self):
        g_eps = self.g_eps
        a,b,h = self.a, self.b, self.h
        no_PML_cells = self.no_PML_cells
        z_measure1 = self.z_measure1
        z_measure2 = self.z_measure2
        self.leq = almost_leq(g_eps)
        self.geq = almost_geq(g_eps)
        self.z_incport_p = close_to_point(self.z_incport, g_eps)        
        self.z_measure1_p = close_to_point(z_measure1, g_eps)
        self.z_measure2_p = close_to_point(z_measure2, g_eps)

        self.left_of_port = lambda r: r[2] < self.z_incport
        self.z_PML1_p = close_to_point(self.z_PML1, g_eps)
        self.z_PML2_p = close_to_point(self.z_PML2, g_eps)
        if not self.short:
            self.sigma_z_fn = PML.fn_1D2fn_3D(
                PML.PML_sigma_fn_1DLR(self.PML_sigma_max,
                                      self.z_PML1, self.z_PML2, self.PML_len), 2)
        else:
            self.z_short_p = close_to_point(self.z_short, g_eps)
            self.sigma_z_fn = PML.fn_1D2fn_3D(
                PML.PML_sigma_fn_1DL(self.PML_sigma_max,
                                     self.z_PML1, self.PML_len), 2)
        self.sigma_fns = Struct(x=lambda r:0, y=lambda r:0, z=self.sigma_z_fn)
        self.zero_p = close_to_point(0, g_eps)
        self.a_p, self.b_p = (close_to_point(xx, g_eps) for xx in (a,b))
        self.z_guide_start_p = close_to_point(self.z_guide_start, g_eps)
        self.z_guide_end_p = close_to_point(self.z_guide_end, g_eps)
        
    def _setup_freefuns(self):
        self.freeE = freeE = cb
        self.freeB = freeE
        self.pml_sys_free = lambda ent: freeE(ent) and self.in_PML_sys(ent)
        self.bonding_sys_free = freeE
        self.bonding_sys_free_vol = self.in_bonding_sys

    def in_bonding_sys(self, ent):
        return not self.in_PML_sys(ent)

    def on_incport(self, ent):
        return N.all(self.z_incport_p(ent.nodeCoords[:,2]))

    def in_PML_sys(self, ent):
        return N.all(self.geq(ent.nodeCoords[:,2], self.z_PML2)
                     | self.leq(ent.nodeCoords[:,2], self.z_PML1))

    def on_z_PML(self, ent):
        return N.all(self.z_PML1_p(ent.nodeCoords[:,2])
                     | self.z_PML2_p(ent.nodeCoords[:,2]))

    def on_PML_bond_surf(self, ent):
        return self.on_z_PML(ent)
        
    def port_free(self, ent):
        ncx = ent.nodeCoords[:,0] ; ncy = ent.nodeCoords[:,1]
        return not (N.all(self.a_p(ncx)) or N.all(self.zero_p(ncx))
                    or N.all(self.b_p(ncy)) or N.all(self.zero_p(ncy)))

    def on_measurement_port1(self, ent):
        return N.all(self.z_measure1_p(ent.nodeCoords[:,2]))

    def on_measurement_port2(self, ent):
        return N.all(self.z_measure2_p(ent.nodeCoords[:,2]))


    def init_mesh(self):
        self.hex_mesh = BrickMesh.Mesh(
            BrickMeshGen.make_rect_cavity_brick_listmesh(
            self.a, self.b, self.c, self.h,
            grid_offset=[0,0,self.mesh_offset]))
        print 'Hex Mesh elements: ', len(self.hex_mesh.elements)

    def _choose_measure1_sys(self):
        return self.pml_system if not self.measure1_on_incport\
               else self.bonding_system

    def _choose_measure2_sys(self):
        return self.pml_system if not self.short \
               else self.bonding_system

    def init_systems(self):
        order = self.order 
        mesh = self.hex_mesh
        print 'oooooooooorder: ', order
        self.bonding_system = self.BondingSystemClass(self.in_bonding_sys)
        self.pml_system = self.PMLBondedSystemClass(
            mesh, BCs=Struct(E=self.pml_sys_free, B=self.pml_sys_free),
            volBCs=Struct(E=self.in_PML_sys, B=self.in_PML_sys),
            disc_orders=self.disc_orders, el_in_disc=self.in_PML_sys)
        self.bonding_system.add_bonded_system(self.pml_system, self.on_PML_bond_surf)
        self.bonding_system.init_discs(
            mesh, BCs=Struct(E=self.bonding_sys_free, B=self.bonding_sys_free),
            volBCs=Struct(E=self.bonding_sys_free_vol, B=self.bonding_sys_free_vol),
            disc_orders=self.disc_orders, el_in_disc=self.in_bonding_sys)
        self.pml_system.set_bond(self.bonding_system.get_bond(self.pml_system))
        self.pml_system.set_sigmas(self.sigma_fns)
        self.systems = (self.bonding_system, self.pml_system)
        self.measure1_sys = self._choose_measure1_sys()
        self.measure2_sys = self._choose_measure2_sys()
        

    def setupSource(self, source_runtime=None, source_no_PML_cells=None):
        if not source_runtime is None: self.source_runtime = source_runtime
        if not source_no_PML_cells is None:
            self.source_no_PML_cells = source_no_PML_cells
        order = self.order
        discs = self.bonding_system.discs
        self.WGB = WGB = WaveGuideTFSFBootstrap(
            self.a,self.b,self.h, order, self.drv_fun,
            no_PML_cells=self.source_no_PML_cells)
        WGB.init_mesh()
        WGB.init_systems()
        WGB.setupSource()
        self.inc_field = inc_field = Feeds.TFSFIncidentField(
            discs, self.on_incport, self.left_of_port)
        inc_field.init_incsurf(self.port_free)
        self.inc_field_info = inc_info = inc_field.get_inc_info()
        self.B_st_dofnos = inc_info.B.dofnos
        self.E_ts_dofnos = inc_info.E.dofnos
        self.C_st = discs.E.matrix.exteriorDerivative(discs.B)[
            N.ix_(self.B_st_dofnos, self.E_ts_dofnos)]
        self.M_mu_st = discs.B.matrix.mass()[
            N.ix_(self.B_st_dofnos, self.B_st_dofnos)]
        self.M_eps_inv_ts = discs.E.matrix.mass()[
            N.ix_(self.E_ts_dofnos, self.E_ts_dofnos)].todia()
        self.M_eps_inv_ts.data[0][:] = 1/self.M_eps_inv_ts.data[0][:]

    def setupLogging(self):
        msys1 = self.measure1_sys
        msys2 = self.measure2_sys
        divisor = self.log_divisor
        self.sm_m2 = sm_m2 = Feeds.BrickSurfaceFieldMatcher()
        self.sm_m1 = sm_m1 = Feeds.BrickSurfaceFieldMatcher()
        self.E_wg = E_wg = Feeds.gen_E_TE01(self.a,self.b)
        sm_m2.initSubdim(msys2.discs.E,
                        self.on_measurement_port2, self.port_free)
        sm_m1.initSubdim(msys1.discs.E,
                         self.on_measurement_port1, self.port_free)
        self.measure_dofArray2 = sm_m2.matchKnown(self.E_wg)
        self.measure_dofArray1 = sm_m1.matchKnown(self.E_wg)
        self.tif2 = tif2 = Feeds.BrickTangentialFuncProjSurfaceIntegral(
            msys2.discs.E, self.on_measurement_port2, self.port_free)
        self.tif1 = tif1 = Feeds.BrickTangentialFuncProjSurfaceIntegral(
            msys1.discs.E, self.on_measurement_port1, self.port_free)
        self.logged_dofnos2 = logged_dofnos2 = tif2.superDOFMap
        self.logged_dofnos1 = logged_dofnos1 = tif1.superDOFMap
        msys2.addLogger('E', dofnos=logged_dofnos2, divisor=divisor)
        msys1.addLogger('E', dofnos=logged_dofnos1, divisor=divisor)
        
    def getResult(self):
        msys1 = self.measure1_sys ; msys2 = self.measure2_sys
        sm_m2 = self.sm_m2 ; sm_m1 = self.sm_m1
        measure_dofArray2 = self.measure_dofArray2
        measure_dofArray1 = self.measure_dofArray1
        nf2 =  N.dot(measure_dofArray2, sm_m2.subDisc.matrix.mass()*(measure_dofArray2))
        nf1 =  N.dot(measure_dofArray1, sm_m1.subDisc.matrix.mass()*(measure_dofArray1))

        ts_modeintg2_n = N.array([
            N.dot(measure_dofArray2, sm_m2.subDisc.matrix.mass()*(dof_vec))
            for dof_vec in msys2.loggedDOFs.E[tuple(self.logged_dofnos2)].vals],
                                N.float64)/nf2
        ts_modeintg1_n = N.array([
            N.dot(measure_dofArray1, sm_m1.subDisc.matrix.mass()*(dof_vec))
            for dof_vec in msys1.loggedDOFs.E[tuple(self.logged_dofnos1)].vals],
                                 N.float64)/nf1
        return Struct(nf2=nf2, nf1=nf1,  ts_modeintg2_n=ts_modeintg2_n,
                      ts_modeintg1_n=ts_modeintg1_n)

    def set_dt(self, dt):
        self.dt = dt
        self.resetHistory()             # note, WGB's dt will be set in resetHistory
        for system in self.systems:
            if self.useLU: system.useLU = True
            print "Setting timestep"
            system.setTimestep(dt)
        print "Done setting timesteps"

    def runSteps(self, n_steps):
        B_steppers = [system.step_B(n_steps) for system in self.systems]
        E_steppers = [system.step_E(n_steps) for system in self.systems]
        C_st = self.C_st
        M_eps_inv_ts = self.M_eps_inv_ts
        M_mu_st = self.M_mu_st
        dt = self.dt
        self.inc_modeintg.append(0.)
        drv_sys = self.bonding_system
        # This is wrong if the system is re-continued
        e_ts = N.zeros_like(drv_sys.dofs.E.dofArray[self.E_ts_dofnos])
        
        for n in xrange(n_steps):
            for bstep in B_steppers: bstep.next()
            drv_sys.dofs.B.dofArray[self.B_st_dofnos] += dt*(C_st*e_ts)

            inc_vals = self.WGB.next_drive()
            e_ts = inc_vals.E_incdofs
            b_st = inc_vals.B_incdofs
            self.inc_modeintg.append(inc_vals.mode_weighted)
            drv_sys.dofs.E.dofArray[self.E_ts_dofnos] += dt*M_eps_inv_ts*(
                C_st.T*(M_mu_st*b_st))

            for estep in E_steppers: estep.next()

    def get_stepsCompleted(self):
        return self.bonding_system.n

    def resetHistory(self):
        self.WGB.set_dt(self.dt)        # Runs WGM.resetHistory() as a side effect
        self.WGB.init_steppers(self.source_runtime)
        self.inc_modeintg = []
        for system in self.systems:
            system.reset_history()
        self.setupLogging()
        for system in self.systems:
            system.log()                # To add an entry for n=0    

    def _calc_z_PML1_goal(self):
        z_measure1 = self.z_measure1
        h = self.h 
        return z_measure1 - h if self.measure1_on_incport else z_measure1

    def _calc_z_PML2_goal(self):
        return self.z_measure2


class WaveGuide2PortBuff(WaveGuide2Port):
    buffspace = 0.

    def _calc_z_PML1_goal(self):
        z_measure1 = self.z_measure1
        h = self.h ; buffspace = self.buffspace
        if buffspace <= 0: return WaveGuide2Port._calc_z_PML1_goal(self)
        else: return z_measure1 - N.max([h, buffspace])

    def _calc_z_PML2_goal(self):
        z_measure2 = self.z_measure2
        buffspace = self.buffspace
        if buffspace <= 0: return WaveGuide2Port._calc_z_PML2_goal(self)
        else: return z_measure2 + buffspace
        
    def _choose_measure1_sys(self):
        buffspace = self.buffspace
        if buffspace <= 0: return WaveGuide2Port._choose_measure1_sys(self)
        else: return self.bonding_system

    def _choose_measure2_sys(self):
        buffspace = self.buffspace
        if buffspace <= 0: return WaveGuide2Port._choose_measure2_sys(self)
        else: return  self.bonding_system

    def set_dims(self, z_incport, z_measure1, z_measure2, buffspace=0, short_after=None):
        """Note z_measure2 > z_incport >= z_measure1 should be true and the
        distance between z_measure1 and z_incport should be zero or a multiple
        of h"""
        assert(buffspace >= 0)
        self.buffspace = buffspace
        WaveGuide2Port.set_dims(
            self, z_incport, z_measure1, z_measure2, short_after=short_after)

class WaveGuide2PortPECels(WaveGuide2Port):
    def set_PECels(self, PEC_elset):
        self.PEC_elset = PEC_elset

    def _setup_freefuns(self):
        self.freeE = freeE = cb
        self.freeB = freeE
        self.pml_sys_free = lambda ent: freeE(ent) and self.in_PML_sys(ent)


    def bonding_sys_free_vol(self, vol):
        return self.in_bonding_sys(vol) and vol.index not in self.PEC_elset

    def bonding_sys_free(self, ent):
        return self.freeE(ent) and not N.any(\
               [con2 in self.PEC_elset for con2 in ent.connect2elem])

class WaveGuide2PortPECelsBuff(WaveGuide2PortBuff):
    def set_PECels(self, PEC_elset):
        self.PEC_elset = PEC_elset

    def _setup_freefuns(self):
        self.freeE = freeE = cb
        self.freeB = freeE
        self.pml_sys_free = lambda ent: freeE(ent) and self.in_PML_sys(ent)


    def bonding_sys_free_vol(self, vol):
        return self.in_bonding_sys(vol) and vol.index not in self.PEC_elset

    def bonding_sys_free(self, ent):
        return self.freeE(ent) and not N.any(\
               [con2 in self.PEC_elset for con2 in ent.connect2elem])
