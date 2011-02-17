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
from NewCode.IOUtils import read_filelog_lines
from analytical_WG_driver import WGT
g_eps = 1e-10                           # Geometrical tollerance

h = 1/8.
a,b = 1., 0.25
dt_div = 4
dt = 1/16/dt_div
order = 3
l_sc = h
z_inc = 0
z_measure = WGT.test_z + z_inc
z_measure_p = close_to_point(z_measure, g_eps)
#runtime = WGT.t_final
runtime = 125
z_PML_goal = z_measure
z_PML1_goal = -1/2
assert (z_PML1_goal <= -l_sc)
no_PML_cells = 10
PML_m = 3
sigma_factor = 1
z_PML = N.ceil(z_PML_goal/h)*h ; PML_len = no_PML_cells*h
z_PML1 = N.floor(z_PML1_goal/h)*h
#c = (z_PML/h + N.abs(z_PML1/h) + no_PML_cells*2)*h
c = (z_PML/h + no_PML_cells + l_sc/h)*h
wave_imp = 1.
PML_sigma_max = 0.8*(PML_m+1)/h/wave_imp*sigma_factor
z_PML2_p = close_to_point(z_PML, g_eps)
z_PML1_p = close_to_point(z_PML1, g_eps)
# z_PML_p = lambda z: z_PML1_p(z) | z_PML2_p(z) 
z_PML_p = lambda z: z_PML2_p(z) 

z_measure1 = z_PML1
z_measure1_p = z_PML1_p


sigma_z_fn2 = lambda r: N.where(
    r.T[2] >= z_PML,
    PML_sigma_max*((N.abs(r.T[2])-z_PML)/PML_len)**PML_m,
    0.)

sigma_z_fn1 = lambda r: N.where(
    r.T[2] <= z_PML1,
    PML_sigma_max*((N.abs(r.T[2] -z_PML1))/PML_len)**PML_m,
    0.)

#sigma_z_fn = lambda r: sigma_z_fn1(r) + sigma_z_fn2(r)

sigma_z_fn = lambda r: sigma_z_fn2(r)


# mesh = BrickMesh.Mesh(
#     BrickMeshGen.make_rect_cavity_brick_listmesh(
#     a,b,c,h,grid_offset=[0,0,z_PML1-PML_len]))
mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(
    a,b,c,h,grid_offset=[0,0,-l_sc]))

print 'Mesh elements: ', len(mesh.elements)

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
constrained = DifferentialForm.allconstrained
z_port = z_inc
z_port_p = close_to_point(z_port, g_eps)
zero_p = close_to_point(0, g_eps)
a_p, b_p = (close_to_point(xx, g_eps) for xx in (a,b))
on_port = lambda ent: N.all(z_port_p(ent.nodeCoords[:,2]))
left_of_port = lambda r: r[2] < z_port
in_PML_sys = lambda ent: N.all(ent.nodeCoords[:,2] >= z_PML) \
             or N.all(ent.nodeCoords[:,2] <= z_PML1)
on_z_PML = lambda ent: N.all(z_PML_p(ent.nodeCoords[:,2])) \
           or N.all(z_PML1_p(ent.nodeCoords[:,2]))
on_PML_bond_surf = on_z_PML                               
in_bonding_sys = lambda ent: not in_PML_sys(ent)

def port_free(ent):
    ncx = ent.nodeCoords[:,0] ; ncy = ent.nodeCoords[:,1]
    return not (N.all(a_p(ncx)) or N.all(zero_p(ncx))
                or N.all(b_p(ncy)) or N.all(zero_p(ncy)))

on_measurement_port = lambda ent: N.all(z_measure_p(ent.nodeCoords[:,2]))
on_measurement_port1 = lambda ent: N.all(z_measure1_p(ent.nodeCoords[:,2]))
freeE = cb
freeB = freeE
pml_sys_free = lambda ent: freeE(ent) and in_PML_sys(ent)
bonding_sys_free = lambda ent: in_bonding_sys(ent) and freeE(ent)

class TestRun(Runners.TestRun):
    drv_fun = staticmethod(lambda *x: 0.)
    useLU = True
    BondingSystemClass = CoupledBondingSystem
    PMLBondedSystemClass = PMLBondedSystem
    sigma_fns = Struct(x=lambda r:0, y=lambda r:0, z=sigma_z_fn)
    
    def __init__(self, mesh, **kwargs):
        self.order = order = kwargs['disc_orders']['E'][0]
        print 'oooooooooorder: ', order
        self.bonding_system = self.BondingSystemClass(in_bonding_sys)
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
        self.pml_system.set_sigmas(self.sigma_fns)
        self.systems = (self.bonding_system, self.pml_system)

    def setupSource(self, input_filename):
        order = self.order
        discs = self.bonding_system.discs
        self.inc_field = inc_field = Feeds.TFSFIncidentField(
            discs, on_port, left_of_port)
        inc_field.init_incsurf(port_free)
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
        self.B_reader = read_filelog_lines(file(input_filename+'_B.numpyraw'))
        self.B_reader.next()            # discard B field at n - 1/2
        self.E_reader = read_filelog_lines(file(input_filename+'_E.numpyraw'))

    def setupLogging(self):
        divisor = self.log_divisor
        from NewCode import PostProc
        self.sm_m = sm_m = Feeds.BrickSurfaceFieldMatcher()
#         self.sm_m1 = sm_m1 = Feeds.BrickSurfaceFieldMatcher()
        self.E_wg = E_wg = Feeds.gen_E_TE01(a,b)
        sm_m.initSubdim(self.pml_system.discs.E, on_measurement_port, port_free)
#         sm_m1.initSubdim(self.pml_system.discs.E, on_measurement_port1, port_free)
        self.measure_dofArray = sm_m.matchKnown(self.E_wg)
#         self.measure_dofArray1 = sm_m1.matchKnown(self.E_wg)
        self.tif = tif = Feeds.BrickTangentialFuncProjSurfaceIntegral(
            self.pml_system.discs.E, on_measurement_port, port_free)
#         self.tif1 = tif1 = Feeds.BrickTangentialFuncProjSurfaceIntegral(
#             self.pml_system.discs.E, on_measurement_port1, port_free)
        self.logged_dofnos = logged_dofnos = tif.superDOFMap
#         self.logged_dofnos1 = logged_dofnos1 = tif1.superDOFMap
        self.pml_system.addLogger('E', dofnos=logged_dofnos, divisor=divisor)
#         self.pml_system.addLogger('E', dofnos=logged_dofnos1, divisor=divisor)
        
    def getResult(self):
        sm_m = self.sm_m
#         sm_m1 = self.sm_m1
        measure_dofArray = self.measure_dofArray
#         measure_dofArray1 = self.measure_dofArray1
        nf =  N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(measure_dofArray))
#         nf1 =  N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(measure_dofArray1))

        ts_modeintg_n = N.array([
            N.dot(measure_dofArray, sm_m.subDisc.matrix.mass().matvec(dof_vec))
            for dof_vec in self.pml_system.loggedDOFs.E[tuple(self.logged_dofnos)].vals],
                                N.float64)/nf        
#         ts_modeintg1_n = N.array([
#             N.dot(measure_dofArray1, sm_m1.subDisc.matrix.mass().matvec(dof_vec))
#             for dof_vec in self.pml_system.loggedDOFs.E[tuple(self.logged_dofnos1)].vals],
#                                  N.float64)/nf1
#         return Struct(nf=nf, ts_modeintg_n=ts_modeintg_n,
#                       ts_modeintg1_n=ts_modeintg1_n)
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
        drv_sys = self.bonding_system
        B_steppers = [system.step_B(n_steps) for system in self.systems]
        E_steppers = [system.step_E(n_steps) for system in self.systems]
        C_st = self.C_st
        M_eps_inv_ts = self.M_eps_inv_ts
        M_mu_st = self.M_mu_st
        
        dt = self.dt
        for n in xrange(n_steps):
            for bstep in B_steppers: bstep.next()
            try:
                e_ts = self.E_reader.next()
                drv_sys.dofs.B.dofArray[self.B_st_dofnos] += dt*C_st.matvec(e_ts)
            except StopIteration: pass

            for estep in E_steppers: estep.next()
            try:
                b_st = self.B_reader.next()
                drv_sys.dofs.E.dofArray[self.E_ts_dofnos] += dt*M_eps_inv_ts.matvec(
                    C_st.T.matvec(M_mu_st.matvec(b_st)))
            except StopIteration: pass

    def get_stepsCompleted(self):
        return self.bonding_system.n

    def resetHistory(self):
        for system in self.systems:
            system.reset_history()
        self.setupLogging()
        for system in self.systems:
            system.log()                # To add an entry for n=0    



input_filebasename =  '+WG_inc'
input_filename = input_filebasename + "%sX%s_h%s_o_%s_dt%s_pml%s" %(
    str(a), str(b), str(1/h), str(order), str(1/dt), str(40))
source_info = pickle.load(file(input_filename+'.pickle'))

E_order = (order, True)
B_order = (order, True)
print 'order: ', order
TRS = TestRun(mesh, disc_orders=dict(E=E_order, B=B_order),
              btype='cohen98')
TRS.log_divisor = dt_div
TRS.setupSource(input_filename)
TRS.set_dt(dt)
no_steps = int(N.ceil(runtime/dt))
TRS.runSteps(no_steps)
resses = {order:{dt_div:TRS.getResult()}}
drv_ts = source_info['ts_modeintg_n'][::dt_div]
resses[order][dt_div]['drv_ts'] = -drv_ts

pickle.dump(resses, file('+brick_4_wg_tfsf_tmp.pickle', 'w'))
f_min, f_max = 0.4375, 2.375
res_ts = resses[order][dt_div].ts_modeintg_n
n_ts = len(res_ts)
#n = 2**int(N.ceil(N.log2(n_ts)))
n = 2**14
assert(n > n_ts)
df = WGT.f_nyq/n*2
res_fs = scipy.fft(res_ts, n=n)
drv_fs = scipy.fft(drv_ts, n=n)
n_f_min, n_f_max = int(N.floor(f_min/df)), int(N.ceil(f_max/df))
tf = (res_fs/drv_fs)[n_f_min:n_f_max]
fr = N.arange(n_f_min, n_f_max)*df
tr = N.arange(0, no_steps+1, dt_div)*dt
