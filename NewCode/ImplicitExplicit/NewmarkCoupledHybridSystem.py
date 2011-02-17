from __future__ import division

from itertools import chain

import numpy as N
from numpy.testing import assert_equal
import scipy

from NewCode.Utilities import Struct
from NewCode.ImplicitExplicit.NewmarkHybridSystem import NewmarkHybridSystem
import NewCode.ImplicitExplicit as IE
import NewCode.ImplicitExplicit.NewmarkHybridMatrices as NHM

def newmark_leapfrog_step(self, no_steps=1):
    C_dc = self.block_matrices.C_dc()
    C_dd = self.block_matrices.C_dd()
    C_ed = self.block_matrices.C_ed()
    C_ee = self.block_matrices.C_ee()
    A_imp = self.merged_matrices.A_imp()
    if self.useLU:
        A_imp.solve = self.merged_matrices.A_imp_LU().solve
    M_dd_diag = self.discs.E.d.matrix.mass().diagonal()
    M_ee_diag = self.discs.E.e.matrix.mass().diagonal()
    M_dd = Struct(matvec=lambda x: M_dd_diag*x,
                   solve=lambda x: x/M_dd_diag)
    M_ee = Struct(matvec=lambda x: M_ee_diag*x,
                   solve=lambda x: x/M_ee_diag)        
    P_dd = self.block_matrices.P_dd()
    P_de = self.block_matrices.P_de()
    P_ee = self.block_matrices.P_ee()
#     A_exp_diag = self.merged_matrices.A_exp().diagonal()
#     A_exp = Struct(matvec=lambda x: A_exp_diag*x,
#                    solve=lambda x: x/A_exp_diag)
    B = self.block_matrices.B_impblocks()
    dt = self.dt
    b = self.implicit_beta
    drv_dofnos, drv_weights, drv_fun = self.drive
    direch = False
    d_o = self.dof_offsets.E
    sa = slice(d_o.a, d_o.b) ; sb = slice(d_o.b, d_o.c) 
    sc = slice(d_o.c, d_o.d) ; sd = slice(d_o.d, d_o.e) 
    slices = Struct(a=sa, b=sb, c=sc, d=sd)
    d_o_B = self.dof_offsets.B
    d_o_E_exp_d = 0; d_o_E_exp_e = d_o.e - d_o.d
    sd_B = slice(d_o_B.d, d_o_B.e) ; se_B = slice(d_o_B.e, self.totalDOFs.B)
    slices_B = Struct(d=sd_B, e=se_B)
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
        else:
            df = self.discs.B[self.direch_group]
            C_p = dd.matrix.exteriorDerivative(df)
            s_direch_B = slices_B[self.direch_group]
        direch = True
        
    for step in xrange(no_steps):
        self.n += 1
        self.newmark_dofArrays[-1][:] = self.newmark_dofArrays[0]
        self.newmark_dofArrays[0][:] = self.newmark_dofArrays[1]
        dm1, d0 = self.newmark_dofArrays[-1], self.newmark_dofArrays[0]
        d0a = d0[sa] ; d0b = d0[sb] ; d0c = d0[sc] ; d0d = d0[sd] 
        y_imp = N.zeros_like(d0[d_o.a:d_o.d])
        y_imp[sa] =                    B.aa.matvec(d0a) + B.ab.matvec(d0b)
        y_imp[sb] = B.ba.matvec(d0a) + B.bb.matvec(d0b) + B.bc.matvec(d0c)
        y_imp[sc] = B.cb.matvec(d0b) + B.cc.matvec(d0c) + B.cd.matvec(d0d)
        y_imp -= A_imp.matvec(dm1[0:d_o.d])

        drvs = [drv_fun(dt, self.n - i) for i in range(3)]
        drv_n = drvs[1]
        y_imp[drv_dofnos] -= drv_weights*(
            b*drvs[0] + (1-2*b)*drvs[1] + b*drvs[2])
        if direch:
            if self.direch_group in self.implicit_groups:
                dp1_p, d0_p, dm1_p = [drv*direch_dofarr for drv in drvs]
                y_imp[s_direch] += B_p.matvec(d0_p) - A_p.matvec(dp1_p + dm1_p)
  
        d0l = self.leapfrog_dofArrays.E
        d0l_B = self.leapfrog_dofArrays.B
        d0l_B[sd_B] -= dt*(C_dc.matvec(d0[sc]) + C_dd.matvec(d0[sd]))
        d0l_B[se_B] -= dt*(C_ed.matvec(d0[sd]) + C_ee.matvec(d0l[d_o_E_exp_e:]))
        if direch:
            if self.direch_group not in self.implicit_groups:
                d0l_B[s_direch_B] -= dt*drv_n*(C_p.matvec(direch_dofarr))
        d0l[d_o_E_exp_d:d_o_E_exp_e] += M_dd.solve(
            P_dd.matvec(d0l_B[sd_B]) + P_de.matvec(d0l_B[se_B]))*dt
        d0l[d_o_E_exp_e:] += M_ee.solve(P_ee.matvec(d0l_B[se_B]))*dt
        self.newmark_dofArrays[1][0:d_o.d] = self.solver.solve_mat_vec(A_imp, y_imp)
        self.log()
        print 'Step %d/%d, drv_fun: %f, total steps: %d, max: %f' % \
                (step+1, no_steps, drvs[0],
                 self.n, N.max(N.abs(self.global_dofarrays.E.all)))



class NewmarkCoupledHybridSystem(NewmarkHybridSystem):
    explicit_groups = ('d', 'e')
    implicit_groups = ('a', 'b', 'c')
    direch_group = 'a'
    
    def init_discs(self, order=1):
        NewmarkHybridSystem.init_discs(self, order=order)
        gff_c, vol_gff_c = self.group_freefuns.c
        gff_d, vol_gff_d = self.group_freefuns.d
        gff_e, vol_gff_e = self.group_freefuns.e
        freeB_d = lambda ent: gff_c(ent) or gff_d(ent)
        self.discs.B = Struct(d=self.DiscretiserModule.setup_PformDiscretiser(
            self.mesh, form=2,order=order,mixed=True,
            freeFun=freeB_d, vol_freeFun=vol_gff_d, btype='cohen98'),
                              e=self.DiscretiserModule.setup_PformDiscretiser(
            self.mesh, form=2,order=order,mixed=True,
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
                k, disc.newDOFs()) for k,disc in discs.iteritems())
            grp_names = sorted(discs.keys())
            dof_lens = [len(bdofs[k].dofArray) for k in grp_names]
            dof_offsets = [0] + N.cumsum(dof_lens[:-1]).tolist()
            tdofs = self.totalDOFs[phys_quant] = sum(dof_lens)
            grp_offsets = self.dof_offsets[phys_quant] = Struct(zip(grp_names, dof_offsets))
            self.global_dofarrays[phys_quant] = g_darr = Struct(all=N.zeros(
                tdofs, dtype))
            for k, start, stop in zip(
                grp_names, dof_offsets, dof_offsets[1:]+[tdofs]):
                bdofs[k].dofArray = g_darr.all[start:stop]
            if phys_quant == 'E':
                try: g_darr.imp = g_darr.all[:grp_offsets.e]
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
