from __future__ import division

import numpy as N
import scipy

from NewCode.Utilities import Struct
from NewCode.ImplicitExplicit.ImplicitExplicit import gen_hybmesh_group_freefuns
from NewCode.ImplicitExplicit.NewmarkCoupledHybridSystem import newmark_leapfrog_step
from NewCode.ImplicitExplicit.HybridMeshNewmarkHybridSystem \
     import HybridMeshNewmarkHybridSystem as HMN

class HybridMeshNewmarkCoupledHybridSystem(HMN):
    def init_group_freefuns(self, on_hbdry, global_freefuns, dead_elements=None):
        HMN.init_group_freefuns(
            self, on_hbdry, global_freefuns.E, dead_elements=dead_elements)
        gff_c, vol_gff_c = self.group_freefuns.c
        gff_d, vol_gff_d = self.group_freefuns.d
        gff_gen_B, tmp = gen_hybmesh_group_freefuns(
            self.meshes, on_hbdry, global_freefuns.B, self.elgroups,
            dead_elements=dead_elements)
        self.B_group_freefuns = Struct(
            d=(lambda ent: gff_c(ent) or gff_d(ent), vol_gff_d),
            e=gff_gen_B.e)
        
    def init_discs(self, order=1, **kwargs):
        HMN.init_discs(self, order, **kwargs)
        gff_d, vol_gff_d = self.B_group_freefuns.d
        gff_e, vol_gff_e = self.B_group_freefuns.e
        
        idm = self.ExplicitDiscretiserModule
        emesh = self.meshes.exp
        self.discs.B = Struct(d=idm.setup_PformDiscretiser(
            emesh, form=2,order=order,mixed=True,
            freeFun=gff_d, vol_freeFun=vol_gff_d, btype='cohen98'),
                              e=idm.setup_PformDiscretiser(
            emesh, form=2,order=order,mixed=True,
            freeFun=gff_e, vol_freeFun=vol_gff_e, btype='cohen98'))
        for disc in self.discs.B.values():
            disc.diagonalise()

    def init_dofs(self, dtype=N.float64):
        self.block_dofs = Struct()
        self.totalDOFs = Struct()
        self.dof_offsets = Struct()
        self.global_dofarrays = Struct()
        for phys_quant, discs in self.discs.iteritems():
            bdofs = self.block_dofs[phys_quant] = Struct((
                k, disc.newDOFs()) for k,disc in discs.iteritems()
                                                         if k != 'c_imp')
            grp_names = sorted(discs.keys())
            try: grp_names.remove('c_imp')
            except ValueError: pass
            dof_lens = [len(bdofs[k].dofArray) for k in grp_names]
            dof_offsets = [0] + N.cumsum(dof_lens[:-1]).tolist()
            tdofs = self.totalDOFs[phys_quant] = sum(dof_lens)
            grp_offsets = self.dof_offsets[phys_quant] = Struct(zip(grp_names, dof_offsets))
            try: grp_offsets.c = grp_offsets.c_exp
            except AttributeError: pass
            self.global_dofarrays[phys_quant] = g_darr = Struct(all=N.zeros(
                tdofs, dtype))
            for k, start, stop in zip(
                grp_names, dof_offsets, dof_offsets[1:]+[tdofs]):
                bdofs[k].dofArray = g_darr.all[start:stop]
            if phys_quant == 'E':
                try: g_darr.imp = g_darr.all[grp_offsets.a:grp_offsets.e]
                except AttributeError: g_darr.imp = g_darr.all
                g_darr.exp = g_darr.all[grp_offsets.d:]

        gdas = self.global_dofarrays
        self.newmark_dofArrays = {-1:gdas.E.imp.copy(),
                                  0:gdas.E.imp.copy(),
                                  1:gdas.E.imp}
        self.leapfrog_dofArrays = Struct(E=gdas.E.exp,
                                         B=gdas.B.all)
        self.solver = self.block_dofs.E.a.solver
        self.reset_history(zero_dofs=False)
        

    def zero(self):
        """ Set all DOFs to 0 """
        for da in self.newmark_dofArrays.itervalues(): da[:] = 0
        self.leapfrog_dofArrays.E[:] = 0
        self.leapfrog_dofArrays.B[:] = 0

        for dof_name,dof in self.block_dofs.E.iteritems():
                assert N.all(dof.dofArray[:] == 0)
    step = newmark_leapfrog_step
