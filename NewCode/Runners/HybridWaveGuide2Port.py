from __future__ import division
import os
from itertools import izip

import numpy as N

from NewCode.Utilities import Struct,  close_to_point
from NewCode.GeomGen.Hybrid import SimpleSurroundedRect, femmesh_tris2ent_nodesets
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Meshes.MeshIO import Femmesh
from NewCode.ImplicitExplicit.HybridMeshNewmarkBondedSystem import \
     HybridMeshNewmarkBondedSystem

from NewCode.Runners.WaveGuide2Port import WaveGuide2PortBuff

class HybridWaveGuide2PortRect(WaveGuide2PortBuff):
    HybridSystemClass = HybridMeshNewmarkBondedSystem
    def init_tetmesh(self, tet_mesh):
        self.tet_mesh = tet_mesh 
        working_dir = tet_mesh.FemmeshDir
        self.bdry_tri_geom = Femmesh.get_femmesh_tris(
            file(working_dir+os.sep+tet_mesh.FemmeshFilename))
        self.tet_material_ent_nodesets = femmesh_tris2ent_nodesets(self.bdry_tri_geom)
        self.bdry_ent_nodes = self.tet_material_ent_nodesets[1]

    def init_hyb_geom(self, tetrange_n, tetrange_p):
        self.hyb_geom = SimpleSurroundedRect(self.hex_mesh)
        self.hyb_geom.set_tet_range(tetrange_n, tetrange_p)
        self.hyb_geom.init_dead_background_element_set()
        self.on_hex_hbdry = self.hyb_geom.on_hbdry_hex

    def on_tet_hbdry(self, ent):
        return tuple(ent.nodes) in self.bdry_ent_nodes \
               and self.freeE(ent)      # Get rid of constrained hbdry tris

    def on_hbdry(self, ent):
        return self.on_tet_hbdry(ent) if ent.meshtype == 'tet' \
               else self.on_hex_hbdry(ent)
        
    def on_hybsys_bond_surf(self, ent):
        ei = self.hybrid_system.elgroups.ei
        con_ei = N.array([con in ei for con in ent.connect2elem], N.bool8)
        return N.any(con_ei) and (not N.all(con_ei)) and (not on_hbdry(ent))

    def in_bonding_sys(self, ent):
        return not self.in_PML_sys(ent) and not \
               N.any([con in self.dead_and_ei_elset for con in ent.connect2elem])

    def in_bonding_sys_vol(self, el):
        return not self.in_PML_sys(el) and not el.index in self.dead_and_ei_elset

    def _setup_freefuns(self):
        def freeE(ent):
            x,y,z = ent.nodeCoords.T
            a_p, b_p = self.a_p, self.b_p,
            zero_p = self.zero_p
            zs_p, ze_p = self.z_guide_start_p, self.z_guide_end_p
            return not (N.all(a_p(x)) or N.all(zero_p(x)) or 
                        N.all(b_p(y)) or N.all(zero_p(y)) or
                        N.all(zs_p(z)) or N.all(ze_p(z)))

        self.freeE = freeE 
        self.freeB = freeE
        self.pml_sys_free = lambda ent: freeE(ent) and self.in_PML_sys(ent)
        self.bonding_sys_free = None    # Must be set using the hybrid system

    def freeE_hyb(self, ent):
        return self.freeE(ent)

    def freeB_hyb(self, ent):
        return self.freeB(ent)
        
    def init_systems(self, reduced_tet_order=0):
        order = self.order 
        hex_mesh = self.hex_mesh
        print 'oooooooooorder: ', order
        self._init_hybsys(reduced_tet_order=reduced_tet_order)
        self._init_pmlsys()
        self._init_bondingsys()
        self.systems = (self.hybrid_system, self.bonding_system, self.pml_system)
        self.measure1_sys = self._choose_measure1_sys()
        self.measure2_sys = self._choose_measure2_sys()

    def _init_bondingsys(self):
        order = self.order 
        mesh = self.hex_mesh
        self.bonding_system = self.BondingSystemClass(self.in_bonding_sys,
                                                      self.in_bonding_sys_vol)
        # Add bonded-systems info to bonding_system
        self.bonding_system.add_bonded_system(self.pml_system, self.on_PML_bond_surf)
        self.bonding_system.add_bonded_system(
            self.hybrid_system, self.on_hybsys_bond_surf)
        self.bonding_sys_free = self.hybrid_system.group_freefuns.e[0]
        self.bonding_system.init_discs(
            mesh, BCs=Struct(E=self.bonding_sys_free, B=self.bonding_sys_free),
            disc_orders=self.disc_orders)#, el_in_disc=self.in_bonding_sys)
        bs_elent = self.bonding_system.discs.E.elements.entity
        bs_elent.dead_elements = self.hyb_geom.dead_background_elements
        # Set bonded-systems bond info
        self.hybrid_system.set_bond(self.bonding_system.get_bond(self.hybrid_system))
        self.pml_system.set_bond(self.bonding_system.get_bond(self.pml_system))

    def _init_pmlsys(self):
        order = self.order 
        mesh = self.hex_mesh
        self.pml_system = self.PMLBondedSystemClass(
            mesh, BCs=Struct(E=self.pml_sys_free, B=self.pml_sys_free),
            volBCs=Struct(E=self.in_PML_sys, B=self.in_PML_sys),
            disc_orders=self.disc_orders, el_in_disc=self.in_PML_sys)
        self.pml_system.set_sigmas(self.sigma_fns)
        
    def _init_hybsys(self,reduced_tet_order):
        order = self.order 
        hex_mesh = self.hex_mesh
        tet_mesh = self.tet_mesh
        self.hybrid_system = self.HybridSystemClass(
            tet_mesh, hex_mesh, implicit_beta=0.25)
        self.hybrid_system.drive= (0,0,lambda *x: 0)
        self.hybrid_system.init_group_freefuns(
            self.on_hbdry, Struct(E=self.freeE_hyb,B=self.freeB_hyb))
        # Bit of a hack...
        self.hybrid_system.elgroups.ei -= self.hyb_geom.dead_background_elements
        self.hybrid_system.init_discs(
            order=order, dead_elements=self.hyb_geom.dead_background_elements,
            reduced_tet_order=reduced_tet_order)
        self.hybrid_system.init_block_matrices()
        self.hybrid_system.init_merged_mats()
        self.hybrid_system.init_dofs()
        self.dead_and_ei_elset = self.hybrid_system.elgroups.ei.union(
            self.hyb_geom.dead_background_elements)

class HybridWaveGuide2PortPECRect(HybridWaveGuide2PortRect):
    def init_hybrid_PEC(self, PEC_label):
        self.hyb_PEC_ent_nodes = self.tet_material_ent_nodesets[PEC_label]
    
    def freeE_hyb(self,ent):
        return self.freeE(ent) and not (tuple(ent.nodes) in self.hyb_PEC_ent_nodes 
                                        if ent.meshtype=='tet' else False)


class HybridWaveGuide2PortPEC(HybridWaveGuide2PortPECRect):
    def init_hyb_geom(self, hyb_geom):
        self.hyb_geom = hyb_geom
        self.hyb_geom.init_dead_background_element_set()
        self.on_hex_hbdry = self.hyb_geom.on_hbdry_hex
