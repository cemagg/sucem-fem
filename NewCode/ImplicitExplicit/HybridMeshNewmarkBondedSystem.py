from __future__ import division

import numpy as N

from NewCode.Utilities import Struct
from NewCode.ImplicitExplicit.HybridMeshNewmarkCoupledHybridSystem \
     import HybridMeshNewmarkCoupledHybridSystem
from NewCode.ImplicitExplicit.HybridMeshNewmarkHybridSystem \
     import HybridMeshNewmarkHybridSystem
from NewCode.ImplicitExplicit.NewmarkBondedSystem \
     import NewmarkBondedSystem


class HybridMeshNewmarkBondedSystem(HybridMeshNewmarkCoupledHybridSystem,
                                    NewmarkBondedSystem):
    explicit_groups = ('c_exp', 'd',)

    def init_discs(self, order=1, **kwargs):
        groupfun_e = self.group_freefuns.e
        del self.group_freefuns.e
        HybridMeshNewmarkHybridSystem.init_discs(self, order=order, **kwargs)
        self.group_freefuns.e = groupfun_e
        gff_c, vol_gff_c = self.group_freefuns.c
        gff_d, vol_gff_d = self.group_freefuns.d
        freeB_d = lambda ent: gff_c(ent) or gff_d(ent)
        self.discs.B = Struct(d=self.ExplicitDiscretiserModule.setup_PformDiscretiser(
            self.meshes.exp, form=2,order=order,mixed=True,
            freeFun=freeB_d, vol_freeFun=vol_gff_d, btype=self.explicit_btype))
        for disc in self.discs.B.values():
            disc.diagonalise()
        
