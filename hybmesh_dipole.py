"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division

from itertools import izip
import numpy as N

from NewCode import Exceptions
from NewCode.Utilities import partial, close_to_point, Struct, \
     in_box, on_box_surf
from NewCode.PostProc import LocatePoints
import NewCode.Mesh as Mesh
from NewCode.GeomGen.Hybrid import SurroundedRect
import NewCode.eMAGUSImport as eMAGUSImport
from NewCode import Waveforms
from NewCode.Meshes.MeshIO import Femmesh
from NewCode import Runners
from NewCode.ImplicitExplicit.HybridMeshNewmarkBondedSystem import \
     HybridMeshNewmarkBondedSystem
from NewCode.DiscretisedSystem.CoupledBonding import \
     CoupledBondingSystem, PMLBondedSystem
from NewCode import DifferentialForm
cb = DifferentialForm.constrained_on_boundary
allfree = DifferentialForm.allfree
g_eps = 1e-10                           # Geometrical tollerance

h = 1/8.
order = 2
dt=1/75.
runtime = 5

#tet_geom_minsize = N.array([0.25,0.25,0.25])
tet_geom_minsize = N.array([0.5,0.5,0.5])
#freespace_minsize = N.array([2.5, 2.5, 2.5])
#tet_geom_minsize = N.array([2,2,2.])
freespace_minsize = N.array([1., 1., 1.])*4
#freespace_minsize = N.array([1, 1., 1.])*h*5



no_PML_cells = 0
sigma_factor = 1
#dipole_pt = N.array([0.,-0.093751,-0.09375])
dipole_pt = N.array([0., -0.03,-0.03])
#dipole_pt = N.array([0.,0.0,0.0])
dipole_dir = N.array([0,0,1.])
test_pts = N.array([[0, 0.5875, 0]])+dipole_pt

geom = SurroundedRect(tet_geom_minsize, freespace_minsize, h)
geom.init_background_mesh(pmlcells=no_PML_cells)
geom.init_dead_background_element_set()
geom.init_geom()

hex_mesh = geom.background_mesh
hex_on_hbdry = geom.on_hbdry_hex

eMAGUSImport.init('workspace')
tet_mesh = Mesh.Mesh(eMAGUSImport.get_listmesh())
print 'Tet_mesh elements: ', len(tet_mesh.elements)
bdry_tri_geom = Femmesh.get_femmesh_tris(
    file('workspace/'+tet_mesh.FemmeshFilename))
bdry_ent_nodes = set([tuple(nds) for nds, mat in
                     izip(bdry_tri_geom.nodes, bdry_tri_geom.material)
                     if mat == 1])
bdry_ent_nodes.update(set([tuple(e_nds) for t_nds in bdry_ent_nodes
                           for e_nds in [ [t_nds[0], t_nds[1]],
                                          [t_nds[0], t_nds[2]],
                                          [t_nds[1], t_nds[2]] ]]))

tet_on_hbdry = lambda ent: tuple(ent.nodes) in bdry_ent_nodes

on_hbdry = lambda ent: tet_on_hbdry(ent) if ent.meshtype == 'tet' \
           else hex_on_hbdry(ent)

PML_len = h*no_PML_cells*N.array([1,1,1], N.float64)
PML_m = 3
wave_imp = 1.
PML_sigma_max = 0.8*(PML_m+1)/h/wave_imp*sigma_factor*N.array([1,1,1], N.float64)
PML_pos = geom.background_mesh_size + geom.offs - h*no_PML_cells
def sigma_u(coord_no, r) :
    u = r.T[coord_no] ; sigma_max = PML_sigma_max[coord_no]
    pos = PML_pos[coord_no] ; plen = PML_len[coord_no]
    return N.where(N.abs(u) >= pos,
                   sigma_max*((N.abs(u)-pos)/plen)**PML_m, 0.)
outside_PML_sysp = in_box(-PML_pos, PML_pos)
on_PML_bdryp = on_box_surf(-PML_pos, PML_pos, g_eps)
on_PML_bdry = lambda ent: N.all([on_PML_bdryp(nc) for nc in ent.nodeCoords])
in_PML_sys = lambda ent: (not N.all([outside_PML_sysp(nc) for nc in ent.nodeCoords])
                          or on_PML_bdry(ent))

sigma_x_fn = partial(sigma_u,0)
sigma_y_fn = partial(sigma_u,1)
sigma_z_fn = partial(sigma_u,2)

a,b,c = geom.background_mesh_size + geom.offs
a_p, b_p, c_p = (close_to_point(xx, g_eps) for xx in (a,b,c))

def freeE(ent):
    x,y,z = ent.nodeCoords.T
    return not (N.all(a_p(x)) or N.all(a_p(-x)) or 
                N.all(b_p(y)) or N.all(b_p(-y)) or
                N.all(c_p(z)) or N.all(c_p(-z)))
freeE = allfree
freeB = freeE
pml_sys_free = lambda ent: in_PML_sys(ent) and freeE(ent)

current_waveform = Waveforms.get_d_gaussian(fc=1, tpr=-60)
drv_fun = current_waveform.D             # RHS proportional to d/dt of current

class TestRun(Runners.TestRun):
    drv_fun = staticmethod(drv_fun)
    useLU = True
    HybridSystemClass = HybridMeshNewmarkBondedSystem
    BondingSystemClass = CoupledBondingSystem
    PMLBondedSystemClass = PMLBondedSystem
    sigma_fns = Struct(x=sigma_x_fn, y=sigma_y_fn, z=sigma_z_fn)
    drive_group = 'a'
    alt_drive_group = 'b'

    def __init__(self, tet_mesh, brick_mesh, BCs, on_hbdry, order, **kwargs):
        self.hybrid_system = self.HybridSystemClass(tet_mesh, brick_mesh, **kwargs)
        self.hybrid_system.init_group_freefuns(on_hbdry, BCs)
        # Bit of a hack...
        self.hybrid_system.elgroups.ei -= geom.dead_background_elements
        self.hybrid_system.init_discs(order=order, dead_elements=geom.dead_background_elements)
        self.hybrid_system.init_block_matrices()
        self.hybrid_system.init_merged_mats()
        self.hybrid_system.init_dofs()
        orders = Struct(E=(order, True), B=(order, True))
        ei = self.hybrid_system.elgroups.ei
        def on_hybsys_bond_surf(ent):
            con_ei = N.array([con in ei for con in ent.connect2elem], N.bool8)
            return N.any(con_ei) and (not N.all(con_ei)) and (not on_hbdry(ent))
        dis_bricks = geom.dead_background_elements
        self.in_bonding_sys = in_bonding_sys = lambda ent:(
            #not in_PML_sys(ent) and 
            not N.any([con in dis_bricks for con in ent.connect2elem]) and
            not N.any([con in ei for con in ent.connect2elem]))
        self.in_bonding_sys_vol = in_bonding_sys_vol = lambda el:(
            #not in_PML_sys(ent) and
            not el.index in ei and
            not el.index in dis_bricks)
#         self.pml_system = self.PMLBondedSystemClass(
#             brick_mesh, BCs=Struct(E=pml_sys_free, B=pml_sys_free),
#             volBCs=Struct(E=in_PML_sys, B=in_PML_sys),
#             disc_orders=orders, el_in_disc=in_PML_sys)
        self.bonding_system = self.BondingSystemClass(in_bonding_sys, in_bonding_sys_vol)
#         self.bonding_system.add_bonded_system(self.pml_system, on_PML_bdry)
        self.bonding_system.add_bonded_system(self.hybrid_system, on_hybsys_bond_surf)
        freeE_e = self.hybrid_system.group_freefuns.e[0]
        bonding_sys_free = freeE_e
        #bonding_sys_free = lambda ent: in_bonding_sys(ent) and cb(ent)
        self.bonding_system.init_discs(
            brick_mesh, BCs=Struct(E=bonding_sys_free, B=bonding_sys_free),
            disc_orders=orders)#, el_in_disc=in_bonding_sys)
        self.bonding_system.discs.E.elements.entity.dead_elements=geom.dead_background_elements
#         self.pml_system.set_bond(self.bonding_system.get_bond(self.pml_system))
        self.hybrid_system.set_bond(self.bonding_system.get_bond(self.hybrid_system))
#         self.pml_system.set_sigmas(self.sigma_fns)
        self.systems = (self.hybrid_system, self.bonding_system, )#self.pml_system)

    def setupSource(self):
        dg = self.drive_group
        dofs = self.hybrid_system.block_dofs.E[dg]
        try: weights, elPerm = dofs.calcProjPointfunRHS_with_elPerm(
            matchfun=lambda r: dipole_dir, r0=dipole_pt)
        except Exceptions.AllZero: weights = N.array([]) ; elPerm = [[], []]
        if len(elPerm[0]) < len(weights) or len(weights) == 0:
            weights_a = weights[elPerm[0]]
            dofs_a = elPerm[1]
            dga = self.alt_drive_group
            dofs = self.hybrid_system.block_dofs.E[dga]
            weights, elPerm = dofs.calcProjPointfunRHS_with_elPerm(
                matchfun=lambda r: dipole_dir, r0=dipole_pt)
            weights_b = weights[elPerm[0]]
            dofs_b = elPerm[1]
            weights = N.hstack((weights_a, weights_b))
            drive_dofnos = N.hstack((dofs_a, dofs_b))
        else:  drive_dofnos = elPerm[1]
        self.hybrid_system.set_driveDOFs(dg, drive_dofnos, weights, self.drv_fun)
        
    def setupLogging(self):
        divisor = self.log_divisor
        hybsys, bondsys = self.hybrid_system, self.bonding_system
        self.log_systems = (hybsys, bondsys)
        #test_elnos, test_el_coords = LocatePoints(hybsys.meshes.imp, test_pts)
        #hybsys.addReconstructedLogger('a', test_elnos, test_el_coords)
        test_elnos, test_el_coords = LocatePoints(bondsys.mesh, test_pts)
        #hybsys.addReconstructedLogger('d', test_elnos, test_el_coords)
        bondsys.addReconstructedLogger('E', test_elnos, test_el_coords)
        assert test_elnos.shape == (1,)
        assert test_elnos[0] not in hybsys.elgroups.ei
        self.recondebugstuff = Struct(test_elnos=test_elnos, test_el_coords=test_el_coords,
                                      test_pts=test_pts)

    def set_dt(self, dt):
        self.dt = dt
        self.resetHistory()
        for system in self.systems:
            if self.useLU: system.useLU = True
            print "Setting timestep"
            system.setTimestep(dt)
        print "Done setting timesteps"

    def getResult(self):
        hybsys, bondsys = self.hybrid_system, self.bonding_system
#        return N.array(hybsys.loggedReconstructed[0].vals)[:,0,] 
#         return (N.array(hybsys.loggedReconstructed[0].vals)[:,0,] +
#                 N.array(bondsys.loggedReconstructed.E[0].vals)[:,0,])
        return N.array(bondsys.loggedReconstructed.E[0].vals)[:,0,]

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


TRS = TestRun(tet_mesh, hex_mesh, BCs=Struct(E=freeE,B=freeB),
              on_hbdry=on_hbdry, implicit_beta=0.25, order=order)
TRS.setupSource()
TRS.log_divisor = 1
TRS.set_dt(dt)
TRS.runSteps(int(N.ceil(runtime/dt))) 

ts = TRS.getResult()

from pylab import plot, xlabel, ylabel, legend, show, title

plot(N.arange(len(ts))*dt, ts[:,0], label='E_x')
plot(N.arange(len(ts))*dt, ts[:,1], label='E_y')
plot(N.arange(len(ts))*dt, ts[:,2], label='E_z')
xlabel('time (s)')
ylabel('E field magnitude at r = [%f, %f, %f]' %
       (test_pts[0][0], test_pts[0][1], test_pts[0][2]) )
legend()

d_edgetable = TRS.hybrid_system.discs.E.d.permuter.globalEntityPermutationTable('edge')
e_edgetable = TRS.bonding_system.discs.E.permuter.globalEntityPermutationTable('edge')
e_facetable = TRS.bonding_system.discs.B.permuter.globalEntityPermutationTable('face')
bond_freeedges = set(i for i,e in enumerate(e_edgetable) if e[0] >=0)
bond_freefaces = set(i for i,e in enumerate(e_facetable) if e[0] >=0)
