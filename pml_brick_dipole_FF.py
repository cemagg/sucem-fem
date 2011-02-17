from __future__ import division

from itertools import izip
import pickle
import numpy as N

from NewCode.Utilities import partial, close_to_point, Struct, \
     in_box, on_box_surf
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode import DifferentialForm, Waveforms, Feeds, Runners
from NewCode.PostProc import LocatePoints
from NewCode.DiscretisedSystem.CoupledBonding import \
     CoupledBondingSystem, PMLBondedSystem, CoupledBondedSystem
from NewCode.Loggers import NTFFCollecter

cb = DifferentialForm.constrained_on_boundary
allfree = DifferentialForm.allfree
allconstrained = DifferentialForm.allconstrained
g_eps = 1e-10                           # Geometrical tollerance

h = 1/8.5
order = 2
dt=1/75
#dt=1/35
runtime = 20
no_PML_cells = 10
dipole_pt = N.array([0,0,0])
dipole_dir = N.array([0,0,1.])


output_filebasename = '+pml_brick_dipole_FF_h' + str(1/h) + '_o' + str(order)\
                      + '_pml' + str(no_PML_cells) + '_dt' + str(1/dt)
drv_waveform = Waveforms.get_d_gaussian(fc=1/2, tpr=-60)

drv_geom_size = N.array([1, 1., 1.])*h*3
drvsys_size = drv_geom_size + 2*h
Gamma_C_size = N.array([2.,2,2])
Gamma_C_corners = (-Gamma_C_size/2, Gamma_C_size/2)
freespace_size = Gamma_C_size + 2*h
sigma_factor = 1

cn = 1                                   # Normalized speed of light

PML_len = h*no_PML_cells*N.array([1,1,1], N.float64)
PML_m = 3
wave_imp = 1.
PML_sigma_max = 0.8*(PML_m+1)/h/wave_imp*sigma_factor*N.array([1,1,1], N.float64)
PML_pos = freespace_size/2

def sigma_u(coord_no, r) :
    u = r.T[coord_no] ; sigma_max = PML_sigma_max[coord_no]
    pos = PML_pos[coord_no] ; plen = PML_len[coord_no]
    return N.where(N.abs(u) >= pos,
                   sigma_max*((N.abs(u)-pos)/plen)**PML_m, 0.)

sigma_x_fn = partial(sigma_u,0)
sigma_y_fn = partial(sigma_u,1)
sigma_z_fn = partial(sigma_u,2)

x_s, y_s, z_s = mesh_size = freespace_size + 2*PML_len
a,b,c = mesh_size - mesh_size/2

brick_listmesh = BrickMeshGen.make_rect_cavity_brick_listmesh(
    x_s, y_s, z_s, h, grid_offset=-mesh_size/2)
brick_mesh = BrickMesh.Mesh(brick_listmesh)
print 'Brick_mesh elements: ', len(brick_mesh.elements)
a_p, b_p, c_p = (close_to_point(xx, g_eps) for xx in (a,b,c))


outside_PML_sysp = in_box(-PML_pos, PML_pos)
on_PML_bdryp = on_box_surf(-PML_pos, PML_pos, g_eps)
on_PML_bdry = lambda ent: N.all([on_PML_bdryp(nc) for nc in ent.nodeCoords])
in_PML_sys = lambda ent: (not N.any([outside_PML_sysp(nc) and not on_PML_bdryp(nc)
                                     for nc in ent.nodeCoords])
                          or on_PML_bdry(ent))

in_drvsysp = in_box(-drvsys_size/2, drvsys_size/2)
on_drvsysp = on_box_surf(-drvsys_size/2, drvsys_size/2, g_eps)
on_drvsys_bdry = lambda ent: N.all([on_drvsysp(nc) for nc in ent.nodeCoords])
in_drvsys = lambda ent: N.all([in_drvsysp(nc) for nc in ent.nodeCoords])


def not_on_mesh_bdry(ent):
    x,y,z = ent.nodeCoords.T
    return not (N.all(a_p(x)) or N.all(a_p(-x)) or 
                N.all(b_p(y)) or N.all(b_p(-y)) or
                N.all(c_p(z)) or N.all(c_p(-z)))

pmlsys_free = lambda ent: in_PML_sys(ent) and not_on_mesh_bdry(ent)
drvsys_freeE = lambda ent: (in_drvsys(ent) or on_drvsys_bdry(ent)) \
               and not_on_mesh_bdry(ent)
drvsys_freeB = drvsys_freeE
drvsys_volfree = in_drvsys
in_bonding_sys = lambda ent: not in_PML_sys(ent) and not in_drvsys(ent)
on_drv_bond_bdry = on_drvsys_bdry
on_pml_bond_bdry = on_PML_bdry


class TestRun(Runners.TestRun):
    useLU = True
    drv_waveform = staticmethod(drv_waveform)
    BondingSystemClass = CoupledBondingSystem
    DrivingSystemClass = CoupledBondedSystem
    PMLBondedSystemClass = PMLBondedSystem
    sigma_fns = Struct(x=sigma_x_fn, y=sigma_y_fn, z=sigma_z_fn)
    def __init__(self, mesh, **kwargs):
        self.order = order = kwargs['disc_orders']['E'][0]
        print 'oooooooooorder: ', order
        self.bonding_system = self.BondingSystemClass(in_bonding_sys)
        self.drv_system = self.DrivingSystemClass(
            mesh, BCs=Struct(E=drvsys_freeE, B=drvsys_freeB),
            # Note volBCs.B must be drvsys_freeE, since drvsys_freeB also
            # returns true if on_scat_bdry(el), which will be true inside the
            # scatterer if it is meshed with one element.
            volBCs=Struct(E=drvsys_freeE, B=drvsys_freeE), 
            disc_orders=kwargs['disc_orders'], el_in_disc=drvsys_volfree)
        self.bonding_system.add_bonded_system(self.drv_system, on_drv_bond_bdry)
        self.pml_system = self.PMLBondedSystemClass(
            mesh, BCs=Struct(E=pmlsys_free, B=pmlsys_free),
            volBCs=Struct(E=in_PML_sys, B=in_PML_sys),
            disc_orders=kwargs['disc_orders'], el_in_disc=in_PML_sys)
        self.bonding_system.add_bonded_system(self.pml_system, on_pml_bond_bdry)
        self.bonding_system.init_discs(
            # Bonding_sytem already applies in_bonding_sys and not on any
            # bonding bdry internally, hence only potential mesh bdry's need to
            # be constrained
            mesh, BCs=Struct(E=not_on_mesh_bdry, B=not_on_mesh_bdry),
            disc_orders=kwargs['disc_orders'], el_in_disc=in_bonding_sys)
        self.pml_system.set_bond(self.bonding_system.get_bond(self.pml_system))
        self.drv_system.set_bond(self.bonding_system.get_bond(self.drv_system))
        self.pml_system.set_sigmas(self.sigma_fns)
        self.systems = (self.bonding_system, self.drv_system, self.pml_system)

    def setupSource(self):
        order = self.order ; drv_sys = self.drv_system
        weights, elPerm = drv_sys.dofs.E.calcProjPointfunRHS_with_elPerm(
            matchfun=lambda r: dipole_dir, r0=dipole_pt)
        drive_dofnos = elPerm[1]
        self.source_info = Struct(weights=weights, elPerm=elPerm,
                                  drive_dofnos=drive_dofnos)
        drv_sys.setDriveDOFs_J(drive_dofnos, weights, self.drv_waveform)
        
    def setupLogging(self):
        print 'Setting up logging'
        E_file = self.E_logfile = file(self.filebasename + '_E.numpyraw', 'wb')
        B_file = self.B_logfile = file(self.filebasename + '_B.numpyraw', 'wb')
        divisor = self.log_divisor
        bondsys = self.bonding_system
        
        self.ntff = NTFFCollecter(bondsys.discs, Gamma_C_corners)
        bondsys.addFileLogger('E', self.ntff.logdofs.E, E_file, divisor)
        bondsys.addFileLogger('B', self.ntff.logdofs.B, B_file, divisor)
#         bondsys.addLogger('E', self.ntff.logdofs.E, divisor)
#         bondsys.addLogger('B', self.ntff.logdofs.B, divisor)
                
    def getResult(self):
        bondsys = self.bonding_system
        self.E_logfile.flush()
        self.B_logfile.flush()
        return Struct(n=bondsys.loggedFileDOFs.E[0].len)
                      
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
    
E_order = (order, True)
B_order = (order, True)
print "Setup TRS"
TRS = TestRun(brick_mesh, disc_orders=dict(E=E_order, B=B_order))
print "Setup Source"
TRS.setupSource()
TRS.filebasename = output_filebasename
print "Setting dt"
TRS.set_dt(dt)
TRS.runSteps(int(N.ceil(runtime/dt)))
res = TRS.getResult()

n_fft = res.n
drv_ts = [drv_waveform(dt, n) for n in xrange(n_fft)]
# import scipy
# E_fft = scipy.fft(res.E_plane.E[:,0])
# n_fft = len(E_fft)
# f_nyq = 1/2/dt
# df = 2*f_nyq/n_fft
# freq = N.arange(n_fft)*df
# drv_fft = scipy.fft(drv_ts)
# E_res = E_fft/drv_fft
# f_min, f_max = 0.4, 2
# fft_min, fft_max = int(N.floor(f_min/df)), int(N.ceil(f_max/df))+1

pickle.dump(Struct(drv_ts=drv_ts, dt=dt, brick_listmesh=brick_listmesh,
                   Gamma_C_corners=Gamma_C_corners, order=order, ntff=Struct(
    log_ents=TRS.ntff.log_ents, globtables=TRS.ntff.globtables)),
            file(output_filebasename +'.pickle', 'w'))

from NewCode.Consts import c0
#plot(freq[fft_min:fft_max]*c0, N.abs(E_res[fft_min:fft_max]))
# temp_discE = DifferentialForm.BrickDiscretiser.setup_PformDiscretiser(
#     brick_mesh, 1, order,freeFun=allfree)



# dofs_1 = tpm.get_direchdofs(1,1)

# direch_dofsys = direch_discE.newDOFs()
# direch_dofsys.dofArray[:] = dofs_1

