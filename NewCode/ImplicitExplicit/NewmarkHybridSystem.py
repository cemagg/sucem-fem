from __future__ import division

import numpy as N
import scipy
from numpy.testing import assert_equal, assert_almost_equal, assert_approx_equal

from NewCode.Utilities import Struct, max_or_scalar
from NewCode.DifferentialForm import BrickDiscretiser, allconstrained, allfree
import NewCode.ImplicitExplicit as IE
import NewCode.ImplicitExplicit.NewmarkHybridMatrices as NHM


class NewmarkHybridSystemBase(object):
    DiscretiserModule = None
    BlockMats = None
    MergedMats = NHM.HybridMergedMats
    group_discretiser_types = Struct(a='imp', b='imp', d='exp', e='exp')

    def initLogs(self):
        self.loggedDOFs = {}
        self.loggedReconstructed = []

    def reset_history(self, zero_dofs=True):
        if zero_dofs: self.zero()
        self.n = 0
        self.initLogs()

    def zero(self):
        """ Set all DOFs to 0 """
        for da in self.dofArrays.itervalues(): da[:] = 0
        for dof_name,dof in self.block_dofs.E.iteritems():
            assert N.all(dof.dofArray[:] == 0)

    def _init_dofs(self, disc_names):
        self.block_dofs = Struct(E=Struct(
            (k,self.discs.E[k].newDOFs()) for k in disc_names))
        dof_lens = [len(self.block_dofs.E[k].dofArray) for k in disc_names]
        dof_offsets = [0] + N.cumsum(dof_lens[:-1]).tolist()
        self.totalDOFs = sum(dof_lens)
        self.dof_offsets  = Struct(zip(disc_names, dof_offsets))
        self.global_dofarrays = g_darr = Struct(all=N.zeros(
            self.totalDOFs, self.block_dofs.E.a.dofArray.dtype))
        for k, start, stop in zip(
            disc_names, dof_offsets, dof_offsets[1:]+[self.totalDOFs]):
            self.block_dofs.E[k].dofArray = g_darr.all[start:stop]
        g_darr.imp = g_darr.all[self.dof_offsets.a:self.dof_offsets.d]
        g_darr.exp = g_darr.all[self.dof_offsets.d:]
        self.dofArrays = {-1:g_darr.all.copy(),
                          0:g_darr.all.copy(),
                          1:g_darr.all}
        self.solver = self.block_dofs.E.a.solver
        self.reset_history(zero_dofs=False)
        
    def init_merged_mats(self):
        self.merged_matrices = self.MergedMats(self.block_matrices)
        
    def set_dt(self, dt):
        self.dt = dt
        self.block_matrices.set_dt(dt)
        try: self.merged_matrices.set_dt(dt)
        except AttributeError: pass

    setTimestep = lambda self, *names, **kwargs: self.set_dt(*names, **kwargs)

    def step(self, no_steps=1):
        A_imp = self.merged_matrices.A_imp()
        if self.useLU:
            A_imp.solve = scipy.sparse.linalg.factorized(A_imp)
        A_exp_diag = self.merged_matrices.A_exp().diagonal()
        A_exp = Struct(matvec=lambda x: A_exp_diag*x,
                       solve=lambda x: x/A_exp_diag)
        B = self.block_matrices.B_allblocks()
        dt = self.dt
        b = self.implicit_beta
        drv_dofnos, drv_weights, drv_fun = self.drive
        direch = False
        d_o = self.dof_offsets
        if hasattr(self, 'direch_disc'):
            dd = self.direch_disc
            df = self.discs.E[self.direch_group]
            M_p = dd.matrix.projectionOnto(df)
            S_p = dd.D().matrix.projectionOnto(df.D())
            if self.group_discretiser_types[self.direch_group] == 'imp':
                A_p = M_p/dt**2 + b*S_p
                B_p = 2*M_p/dt**2 - (1-2*b)*S_p
            else:
                A_p = M_p/dt**2
                B_p = 2*M_p/dt**2 - S_p
            direch_dofarr = self.direch_dofs.dofArray
            no_a = A_p.shape[0]
            direch = True

        for step in xrange(no_steps):
            self.n += 1
            self.dofArrays[-1][:] = self.dofArrays[0]
            self.dofArrays[0][:] = self.dofArrays[1]
            dm1, d0 = self.dofArrays[-1], self.dofArrays[0]
            sa = slice(d_o.a, d_o.b) ; sb = slice(d_o.b, d_o.c) 
            sc = slice(d_o.c, d_o.d) ; sd = slice(d_o.d, d_o.e) 
            se = slice(d_o.e, self.totalDOFs) ; 
            s_imp = slice(d_o.a,d_o.d) ; s_exp = slice(d_o.d,self.totalDOFs)
            d0a = d0[sa] ; d0b = d0[sb] ; d0c = d0[sc] ; d0d = d0[sd] ; d0e = d0[se]
            y = N.zeros_like(d0)
            y[sa] =                    B.aa*(d0a) + B.ab*(d0b)
            y[sb] = B.ba*(d0a) + B.bb*(d0b) + B.bc*(d0c)
            y[sc] = B.cb*(d0b) + B.cc*(d0c) + B.cd*(d0d)
            y[sd] = B.dc*(d0c) + B.dd*(d0d) + B.de*(d0e)
            y[se] = B.ed*(d0d) + B.ee*(d0e)
            dm1_imp = dm1[s_imp] ; dm1_exp = dm1[s_exp]
            y_imp = y[s_imp] ; y_exp = y[s_exp]
            y_imp -= A_imp*(dm1_imp)
            y_exp -= A_exp*(dm1_exp)

            drvs = [drv_fun(dt, self.n - i) for i in range(3)]
            y[drv_dofnos] -= drv_weights*(
                b*drvs[0] + (1-2*b)*drvs[1] + b*drvs[2])
            if direch:
                dp1_p, d0_p, dm1_p = [drv*direch_dofarr for drv in drvs]
                slices = dict(a=sa,b=sb,c=sc,d=sd,e=se)
                y[slices[self.direch_group]] += B_p*(d0_p) \
                                                - A_p*(dp1_p + dm1_p)
            print 'Step %d/%d, drv_fun: %f, total steps: %d, max: %f' % \
                    (step+1, no_steps, drvs[0],
                     self.n, N.max(N.abs(self.dofArrays[1])))
            self.dofArrays[1][s_imp] = self.solver.solve_mat_vec(A_imp, y_imp)
            self.dofArrays[1][s_exp] = self.solver.solve_mat_vec(A_exp, y_exp)
            self.log()

    def log(self):
        for (group, dofNos), log in self.loggedDOFs.items():
            if self.n % log.div == 0:
                log.vals.append(self.block_dofs.E[group].dofArray[list(dofNos)].copy())
        for rLog in self.loggedReconstructed:
            if self.n % rLog.div == 0:
                dofs = self.block_dofs.E[rLog.group]
                rLog.vals.append(dofs.recon_fromparms(
                        rLog.reconparms, dofs.dofArray))

    def addLogger(self, group, dofnos, divisor=1):
        """
        Log values of DOF # dofnos values for discretiser discName every divisor timesteps
        """
        assert max(dofnos) < self.discs.E[group].totalDOFs
        self.loggedDOFs[group, tuple(dofnos)] = Struct(div=divisor, vals=[])

    def addReconstructedLogger(self, group, log_elnos, log_el_coords, divisor=1):
        self.loggedReconstructed.append(Struct(
            reconparms=list(self.block_dofs.E[group].calc_reconparms(
            log_elnos, log_el_coords)), div=divisor, vals=[], group=group))


class NewmarkHybridSystem(NewmarkHybridSystemBase):
    DiscretiserModule = BrickDiscretiser
    BlockMats = NHM.HybridBlockMats
    direch_group = 'a'

    def __init__(self, mesh, implicit_beta=0.25):
        self.mesh = mesh
        self.implicit_beta = implicit_beta
        self.n = 0
        
    def init_elgroups(self, imp_elset):
        self.elgroups = IE.gen_elgroups(self.mesh, imp_elset)
        self.imp_elset = imp_elset

    def init_group_freefuns(self, global_freefun):
        self.global_freefun = global_freefun
        self.group_freefuns = IE.gen_group_freefuns(
            self.mesh, self.elgroups, global_freefun)

    def init_discs(self, order=1):
        self.order = order
        self.discs = Struct(E=Struct((
            k, self.DiscretiserModule.setup_PformDiscretiser(
            self.mesh, form=1,order=order,mixed=True,freeFun=ff,
            vol_freeFun=ff_vol,  btype='cohen98'))
                                     for k, (ff, ff_vol)
                                     in self.group_freefuns.iteritems()))
        for disc in self.discs.E.values():
            disc.diagonalise()

    def init_dofs(self):
        disc_names = sorted(self.discs.E.keys())
        self._init_dofs(disc_names)
        
    def init_block_matrices(self):
        self.block_matrices = self.BlockMats(
            self.discs, self.implicit_beta, self.imp_elset)

    def set_driveDOFs(self, dofnos, weights, drv_fun):
        """Set global dof numbers, weights and drive function for {f}

        Idea is to have per dof-group speficiation or similar in the
        future. Also, for now the correct value of beta is used only if the
        drive DOFs are in group a or b.
        """
        assert(N.isscalar(weights) or len(weights) == len(dofnos))
        assert(max_or_scalar(dofnos) < self.totalDOFs)
        self.drive = (dofnos, weights, drv_fun)

    def set_direchBCs(self, DirechBCs, DirechVolBCs=allconstrained):
        order = self.order
        self.direch_disc = self.DiscretiserModule.setup_PformDiscretiser(
            self.mesh, form=1, order=order, mixed=True, btype='cohen98',
            freeFun=DirechBCs, vol_freeFun=DirechVolBCs)
        self.direch_dofs = self.direch_disc.newDOFs()
        self.direch_disc.diagonalise()



