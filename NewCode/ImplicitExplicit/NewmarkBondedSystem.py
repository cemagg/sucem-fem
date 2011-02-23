from __future__ import division

import numpy as N

from NewCode.Utilities import Struct
from NewCode.ImplicitExplicit.NewmarkCoupledHybridSystem import NewmarkCoupledHybridSystem
from NewCode.ImplicitExplicit.NewmarkHybridSystem import NewmarkHybridSystem

class NewmarkBondedSystem(NewmarkCoupledHybridSystem):
    explicit_groups = ('d',)
    writestuff_exit = False
    def init_discs(self, order=1):
        groupfun_e = self.group_freefuns.e
        del self.group_freefuns.e
        NewmarkHybridSystem.init_discs(self, order=order)
        self.group_freefuns.e = groupfun_e
        gff_c, vol_gff_c = self.group_freefuns.c
        gff_d, vol_gff_d = self.group_freefuns.d
        freeB_d = lambda ent: gff_c(ent) or gff_d(ent)
        self.discs.B = Struct(d=self.DiscretiserModule.setup_PformDiscretiser(
            self.mesh, form=2,order=order,mixed=True,
            freeFun=freeB_d, vol_freeFun=vol_gff_d, btype='cohen98'))
        for disc in self.discs.B.values():
            disc.diagonalise()

    def set_bond(self, bond):
        self.bond = bond

    def get_bond_data(self):
        return Struct(E_disc=self.discs.E.d, E_dofs=self.block_dofs.E.d)

    def step_B(self, no_steps=1):
        C_dc = self.block_matrices.C_dc()
        C_dd = self.block_matrices.C_dd()
        direch = False
        dt = self.dt        
        d_o = self.dof_offsets.E
        sa = slice(d_o.a, d_o.b) ; sb = slice(d_o.b, d_o.c) 
        sc = slice(d_o.c, d_o.d) ; sd = slice(d_o.d, self.totalDOFs.E) 
        if hasattr(self, 'direch_disc'):
            dd = self.direch_disc
            if self.direch_group in self.explicit_groups:
                df = self.discs.B[self.direch_group]
                C_p = dd.matrix.exteriorDerivative(df)
                direch = True
        drv_dofnos, drv_weights, drv_fun = self.drive

        for step in xrange(no_steps):
            self.n += 1
            self.newmark_dofArrays[-1][:] = self.newmark_dofArrays[0]
            self.newmark_dofArrays[0][:] = self.newmark_dofArrays[1]
            dm1, d0 = self.newmark_dofArrays[-1], self.newmark_dofArrays[0]
            drv_n = drv_fun(dt, self.n-1)
            d0l_B = self.leapfrog_dofArrays.B
            d0l_B -= dt*(C_dc*(d0[sc]) + C_dd*(d0[sd]))
            if direch:
                if self.direch_group in self.explicit_groups:
                    d0l_B -= dt*drv_n*(C_p*(direch_dofarr))
            yield

    def step_E(self, no_steps=1):
        C_dd = self.block_matrices.C_dd()
        M_dd_diag = self.discs.E.d.matrix.mass().diagonal()        
        M_dd = Struct(matvec=lambda x: M_dd_diag*x,
                       solve=lambda x: x/M_dd_diag)
        Mb_dd_diag = self.discs.B.d.matrix.mass().diagonal()
        Mb_dd = Struct(matvec=lambda x: Mb_dd_diag*x,
                       solve=lambda x: x/Mb_dd_diag)
        A_imp = self.merged_matrices.A_imp()
        if self.useLU:
            A_imp.solve = self.merged_matrices.A_imp_LU().solve
        B = self.block_matrices.B_impblocks()
        dt = self.dt
        b = self.implicit_beta
        drv_dofnos, drv_weights, drv_fun = self.drive
        direch = False
        d_o = self.dof_offsets.E
        sa = slice(d_o.a, d_o.b) ; sb = slice(d_o.b, d_o.c) 
        sc = slice(d_o.c, d_o.d) ; sd = slice(d_o.d, self.totalDOFs.E) 
        slices = Struct(a=sa, b=sb, c=sc, d=sd)
        if hasattr(self, 'direch_disc'):
            dd = self.direch_disc
            direch_dofarr = self.direch_dofs.dofArray
            if self.direch_group in self.implicit_groups:
                df = self.discs.E[self.direch_group]
                M_p = dd.matrix.projectionOnto(df)
                S_p = dd.D().matrix.projectionOnto(df.D())
                A_p = M_p/dt**2 + b*S_p
                B_p = 2*M_p/dt**2 - (1-2*b)*S_p
                no_a = A_p.shape[0]
                s_direch = slices[self.direch_group]
                direch = True

        for step in xrange(no_steps):
            dm1, d0 = self.newmark_dofArrays[-1], self.newmark_dofArrays[0]
            d0a = d0[sa] ; d0b = d0[sb] ; d0c = d0[sc] ; d0d = d0[sd] 
            y_imp = N.zeros_like(d0[d_o.a:d_o.d])
            y_imp[sa] =                    B.aa*(d0a) + B.ab*(d0b)
            y_imp[sb] = B.ba*(d0a) + B.bb*(d0b) + B.bc*(d0c)
            y_imp[sc] = B.cb*(d0b) + B.cc*(d0c) + B.cd*(d0d)
            y_imp -= A_imp*(dm1[0:d_o.d])

            drvs = [drv_fun(dt, self.n - i) for i in range(3)]
            drv_n = drvs[1]
            y_imp[drv_dofnos] -= drv_weights*(
                b*drvs[0] + (1-2*b)*drvs[1] + b*drvs[2])
            if direch:
                if self.direch_group in self.implicit_groups:
                    dp1_p, d0_p, dm1_p = [drv*direch_dofarr for drv in drvs]
                    y_imp[s_direch] += B_p*(d0_p) - A_p*(dp1_p + dm1_p)
            self.newmark_dofArrays[1][0:d_o.d] = A_imp.solve(y_imp)
            if hasattr(self, 'writestuff_n'):
                if self.n == self.writestuff_n:
                    from NewCode.IOUtils import mm_write
                    import pickle
                    print "Writing A_imp"
                    mm_write(file(self.writefile+'.mtx', 'w'), A_imp)
                    pickle.dump(Struct(dofarrays=self.newmark_dofArrays,
                                       imp_len=d_o.d, RHS=y_imp),
                                file(self.writefile+'-arrays.pickle', 'w'))
                    if self.writestuff_exit:
                        import sys
                        sys.exit()
            #self.solver.solve_mat_vec(A_imp, y_imp)

            C_bm = self.bond.C_bm ; Mb_b = self.bond.Mb_b
            b_b = self.bond.b_b
            d0l = self.leapfrog_dofArrays.E
            d0l_B = self.leapfrog_dofArrays.B
            d0l += dt*M_dd.solve(
                C_dd.T*(Mb_dd*(d0l_B)) + C_bm.T*(Mb_b*(b_b)))

            self.log()
            # print 'Step Newmark Bonded %d/%d, drv_fun: %f, total steps: %d, max: %f' % \
#                   (step+1, no_steps, drvs[0],
#                    self.n, N.max(N.abs(self.global_dofarrays.E.all)))
            yield
