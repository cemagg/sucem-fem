from __future__ import division

import numpy as N
import scipy

from NewCode.Utilities import Struct, max_or_scalar
from NewCode.DifferentialForm import BrickDiscretiser, Discretiser, allconstrained, allfree
from NewCode.Integration import TriIntegrator, TetProdIntegrator

from NewCode.ImplicitExplicit.NewmarkHybridSystem import NewmarkHybridSystemBase
from NewCode.ImplicitExplicit.ImplicitExplicit import gen_hybmesh_group_freefuns
from NewCode.ImplicitExplicit.HybridMeshNewmarkHybridMatrices import \
     HybridMeshHybridBlockMats

class HybridMeshNewmarkHybridSystem(NewmarkHybridSystemBase):
    ExplicitDiscretiserModule = BrickDiscretiser
    ImplicitDiscretiserModule = Discretiser
    implicit_btype = None
    explicit_btype = 'cohen98'
    BlockMats = HybridMeshHybridBlockMats
    implicit_groups = ('a','b', 'c_imp')
    explicit_groups = ('c_exp', 'd', 'e')
    
    def __init__(self, imp_mesh, exp_mesh, implicit_beta=0.25):
        self.meshes = Struct(imp=imp_mesh, exp=exp_mesh)
        self.implicit_beta = implicit_beta
        self.n = 0
    
        
    def init_group_freefuns(self, on_hbdry, global_freefun, dead_elements=None):
        """Initialise freefuns for the DOF groups

        Input Params
        ============
        global_freefun -- Freedom function for whole problem
        on_hbdry -- Function; returns True for geom entities on hybrid boundary


        Initialises Attributes:
        =======================
        global_freefun -- Same as input
        group_freefuns -- Freedom functions for groups a-e
        on_hbdry -- True for geom entities on hybrid boundary
        """
        self.global_freefun = global_freefun
        self.group_freefuns, self.elgroups = gen_hybmesh_group_freefuns(
            self.meshes, on_hbdry, global_freefun, dead_elements=dead_elements)
        self.on_hbdry = on_hbdry

    def _set_exp_intg(self, disc, order):
        disc.diagonalise()

    def _set_imp_intg(self, disc, order):
        try:
            disc.setIntegrationRule(order*2)
            disc.D().setIntegrationRule(order*2)
        except KeyError:
            disc.set_integrator(TetProdIntegrator)
            disc.D().set_integrator(TetProdIntegrator)
            disc.setIntegrationRule(order*2)
            disc.D().setIntegrationRule(order*2)
        disc.setFaceIntegrationRule(TriIntegrator(order*2))    

    def init_discs(self, order=1, dead_elements=None, reduced_tet_order=0):
        """Initialise discretisers for the various DOF groupings
        """
        
        self.order = order
        self.mtet_order = mtet_order = order*2-1 # mtet for Matching Tet
        self.impl_order = iorder = max(order - reduced_tet_order, 1)
        self.discs = Struct(E=Struct())
        idm = self.ImplicitDiscretiserModule ; edm = self.ExplicitDiscretiserModule
        imsh = self.meshes.imp ; emsh = self.meshes.exp
        gff = self.group_freefuns
        self.discs.E.a = idm.setup_PformDiscretiser(
            imsh, form=1, order=iorder, mixed=True, freeFun=gff.a[0],
            vol_freeFun=gff.a[1])
        self.discs.E.b = idm.setup_PformDiscretiser(
            imsh, form=1, order=iorder, mixed=True, freeFun=gff.b[0],
            vol_freeFun=gff.b[1])
        self.discs.E.c_imp = idm.setup_PformDiscretiser(
            imsh, form=1, order=mtet_order, mixed=False, freeFun=gff.c[0],
            vol_freeFun=gff.c[1])
        self.discs.E.c_exp = edm.setup_PformDiscretiser(
            emsh, form=1, order=order, mixed=True, freeFun=gff.c[0],
            vol_freeFun=gff.c[1], btype=self.explicit_btype, dead_elements=dead_elements)
        self.discs.E.d = edm.setup_PformDiscretiser(
            emsh, form=1, order=order, mixed=True, freeFun=gff.d[0],
            vol_freeFun=gff.d[1], btype=self.explicit_btype, dead_elements=dead_elements)
        # NewmarkBondedSystem based hybrids do not need the "e" group discretiser
        try: self.discs.E.e = edm.setup_PformDiscretiser(
            emsh, form=1, order=order, mixed=True, freeFun=gff.e[0],
            vol_freeFun=gff.e[1], btype=self.explicit_btype, dead_elements=dead_elements)
        except AttributeError: pass
        
        for k in self.implicit_groups:
            self._set_imp_intg(self.discs.E[k], mtet_order)
        for k in self.explicit_groups:
            self._set_exp_intg(self.discs.E[k], order)

    def init_dofs(self):
        disc_names = set(self.implicit_groups + self.explicit_groups) \
                     - set(['c_imp'])
        self._init_dofs(sorted(disc_names))
        self.dof_offsets.c = self.dof_offsets.c_exp
        
    def init_block_matrices(self):
        self.block_matrices = self.BlockMats(
            self.discs, self.implicit_beta, self.on_hbdry)

    def set_direchBCs(self, DirechBCs, DirechVolBCs=allconstrained):
        order = self.order ; mtet_order = self.mtet_order
        group = self.direch_group
        gtype = self.group_discretiser_types[group]
        if gtype == 'exp':
            self.direch_disc = self.ExplicitDiscretiserModule.setup_PformDiscretiser(
                self.meshes.exp, form=1, order=order, mixed=True, btype='cohen98',
                freeFun=DirechBCs, vol_freeFun=DirechVolBCs)
            self._set_exp_intg(self.direch_disc, order)
        elif gtype == 'imp':
            self.direch_disc = self.ImplicitDiscretiserModule.setup_PformDiscretiser(
                self.meshes.imp, form=1, order=order, mixed=True,
                freeFun=DirechBCs, vol_freeFun=DirechVolBCs)
            self._set_imp_intg(self.direch_disc, mtet_order)
        else: raise ValueError(
            'Group type must be "imp" for implicit or "exp" for explicit')
        self.direch_dofs = self.direch_disc.newDOFs()

    def set_driveDOFs(self, group, dofnos, weights, drv_fun):
        """Set global dof numbers, weights and drive function for {f}
        """
        self.direch_group = group
        assert(N.isscalar(weights) or len(weights) == len(dofnos))
        try: dofnos += self.dof_offsets.E[group]
        except AttributeError: dofnos += self.dof_offsets[group]
        assert(max_or_scalar(dofnos) < self.totalDOFs)
        self.drive = (dofnos, weights, drv_fun)

