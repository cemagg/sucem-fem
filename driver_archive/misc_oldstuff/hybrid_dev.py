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
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Utilities import Struct, partial, close_to_point, ArrayMemory
from NewCode import SubDimMesh, DifferentialForm, SystemMatrix, PostProc
from NewCode.DifferentialForm import Discretiser, SubDimDiscretiser,SubDimDiscretiserEntities
from NewCode.DifferentialForm import DiscretiserDOFs
from NewCode.DiscretisedSystem import CoupledFirstOrderSystemBDirechlet
from NewCode.DiscretisedSystem import CurlCurlNewmark
from NewCode.DiscretisedSystem import CurlCurlNewmarkDirechletCoupled
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets

from wg_dof_matching_stuff import RoughWGMatcher

mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh=mesh4

print 'Mesh elements: ', len(mesh.elements)

Discretiser.BasePformDiscretiser.defaultIntegrationOrder = 4

class NewmarkCoupledDirechletDriver(object):
    def __init__(self, newmark_sys, coupled_sys, newmarkDirechBC):
        assert (coupled_sys.mesh is newmark_sys.mesh)
        coupled_sys.setNewmarkCoupledSys(newmark_sys)
        self.dt = coupled_sys.dt
        self.newmarkSys = newmark_sys
        self.coupledSys = coupled_sys
        self.direchSys = CurlCurlNewmark(newmark_sys.mesh, order=newmark_sys.order,
                                         BC=newmarkDirechBC, useQ=newmark_sys.hasQ)
        self.direchSys.dof_mem = ArrayMemory(self.direchSys.dofs.dofArray, mem_len=3)
        self.newmarkSys.setDirechSystem(self.direchSys)
        self.coupledDofMap = self.createDofMap(self.direchSys, coupled_sys)
        self.newmarkSys.setTimestep(self.dt)

    def step(self, no_steps):
        for i in range(no_steps):
            print "Hybrid Timestep %d of %d" % (i+1, no_steps)
            self.coupledSys.step(1)
            self.direchSys.dof_mem.push(
                self.coupledSys.dofs.E.dofArray[self.coupledDofMap])
            self.newmarkSys.step(1)
    def getDirechletDofs(self):
        return self.dof_mem

    @staticmethod
    def createDofMap(direch_sys, coupled_sys):
        dofmap = {}
        disc_d, disc_c = direch_sys.disc, coupled_sys.discs.E
        p_d, p_c = (disc_d.permuter.permuteElementEntities,
                    disc_c.permuter.permuteElementEntities)
        for el_d, el_c in izip(disc_d.elements, disc_c.elements):
            assert (N.all(el_d.nodes == el_c.nodes))
            d_perm, c_perm = (p_d(el_d, remove_constrained=False),
                              p_c(el_c, remove_constrained=False))
            for ent_name in d_perm:
                dp, cp = d_perm[ent_name], c_perm[ent_name]
                assert (dp.shape == cp.shape)
                for dofno_d, dofno_c in izip(dp.flat, cp.flat):
                    if dofno_d < 0: continue
                    assert (dofno_c > -1)
                    dofmap[dofno_d] = dofno_c

        coupled_dof_map = N.zeros(len(dofmap), N.int32)
        for dofno_d, dofno_c in dofmap.iteritems():
            coupled_dof_map[dofno_d] = dofno_c
        assert (coupled_dof_map.shape == direch_sys.dofs.dofArray.shape)
        return coupled_dof_map


cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

wgm = RoughWGMatcher(mesh)
eps = wgm.eps
onPort = wgm.onPort
a,b = wgm.a, wgm.b
hybrid_boundary_z = 5.
explicit_direch_freeE = lambda edge: edge.index in wgm.work_edgeno_set

# On the x or y walls of the waveguide
on_wg_wall = lambda ent: not N.all([
    wgm.edge1_x(x) or wgm.edge2_x(x) or wgm.edge3_y(y) or wgm.edge4_y(y)
    for (x,y,z) in edge.nodeCoords])
on_hybrid_boundary_z = close_to_point(hybrid_boundary_z, eps)
# On the boundary defined between the two regions
on_hybrid_boundary = lambda ent: N.all(on_hybrid_boundary_z(ent.nodeCoords[:,2]))
# In the implicit region, including the wg walls and the hybrid boundary
in_implicit_region_z = lambda z: N.logical_or(z > hybrid_boundary_z,
                                               on_hybrid_boundary_z(z))
in_implicit_region = lambda ent: N.all(in_implicit_region_z(ent.nodeCoords[:,2]))
in_explicit_region_z = lambda z: N.logical_or(z < hybrid_boundary_z,
                                              on_hybrid_boundary_z(z))
in_explicit_region = lambda ent: N.all(in_explicit_region_z(ent.nodeCoords[:,2]))

# All E functions in explicit region, excluding the imp/exp boundary and wg walls
explicit_freeE = lambda ent: (cb(ent)
                              and in_explicit_region(ent)
                              and (not on_hybrid_boundary(ent)))
# All B functions in explicit, inc the imp/exp bdry and wg port, excl wg walls
explicit_freeB = lambda ent: (explicit_freeE(ent)
                             or onPort(ent)
                             or on_hybrid_boundary(ent))
# Implicit E funcs, including the imp/exp bdry, excluding wg walls
implicit_freeE = lambda ent: cb(ent) and in_implicit_region(ent)

# Impl E funcs that update hybrd bdry implicit B fns
implicitE_updateding_explicitB = lambda ent: (
    on_hybrid_boundary(ent) and implicit_freeE(ent))

# Expl E funcs that act as inhomog direchlet for Implicit E

hybrid_bdry_face_nos = N.array([i for i, edge in enumerate(mesh.faces)
                             if on_hybrid_boundary(edge)], N.int32)
hybrid_bdry_element_nos = mesh.faces[hybrid_bdry_face_nos].connect2elem.flatten()
explicit_bdry_element_nos = [i for i in hybrid_bdry_element_nos
                             if in_explicit_region_z(mesh.elements[i].midpoint[2])]
implicit_direch_edgenos = N.array([i for i in N.unique(
    mesh.elements[explicit_bdry_element_nos].edgenos.flat)
                          if not on_hybrid_boundary(mesh.edges[i])], N.int32)
implicit_direch_free = lambda edge: (edge.index in implicit_direch_edgenos
                                     and cb(edge))
dt = 0.01
explicit_orders = {'E':(1,True), 'B':(1,True)}

explicit_system = CoupledFirstOrderSystemBDirechlet(
    mesh, dt=dt, BCs=Struct(E=explicit_freeE, B=explicit_freeB),
    disc_orders=explicit_orders)
explicit_system.setDirechBCs(Struct(E=explicit_direch_freeE, B=lambda x: False))
explicit_system.direchSys.dofs.E.dofArray[:] = wgm.dofs

T, omega = 0.5, 5.523599
def drv_fun(dt, n):
    t = dt*n
    return (1. - N.exp(-(t/2/T)**2))*N.sin(omega*t)
explicit_system.drv_fun = drv_fun

implicit_system = CurlCurlNewmarkDirechletCoupled(mesh, order=1, BC=implicit_freeE,
                                                  useQ=False)

hybrid_system = NewmarkCoupledDirechletDriver(
    implicit_system, explicit_system, newmarkDirechBC=implicit_direch_free)

#hybrid_system.step(1000)
