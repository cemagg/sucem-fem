from itertools import izip
import numpy as N
import scipy
from NewCode.MatrixUtils import MatrixSolver
from NewCode import PostProc

class PformDiscretiserBase(object):
    useLU = False
    def __init__(self, discretiser, dtype):
        self.disc = discretiser
        self.dofArray = discretiser.newDOFsArray(dtype)
        self.dtype = dtype
        self.matrix = discretiser.matrix
        self.solver = MatrixSolver()
        
    def zero(self):
        self.dofArray[:] = 0

    def reconstruct(self, elnos, local_coords):
        return self.recon_fromparms(
            self.calc_reconparms(elnos, local_coords), self.dofArray)

    @staticmethod
    def recon_fromparms(reconparms, dofs):
        return N.array([
            (pvals.T*dofs[perm_g]).sum(axis=1)
            for perm_g, pvals in reconparms
            ], dofs.dtype)

    def calc_reconparms(self, elnos, local_coords):
        elements = self.disc.elements
        for el, coord in izip(elements[elnos], local_coords):
            assert (el.noLocalCoords, ) == coord.shape
            pl, pg = el.permutation()
            yield pg, el.physValsAtPoints([coord])[pl,0]

    def _solveForRHS(self, RHS):
        M = self.matrix.mass_LU() if self.useLU else self.matrix.mass()
        self.dofArray[:] = self.solver.solve_mat_vec(M, RHS)
        
    def matchFunction(self, matchfun):
        self._solveForRHS(self.calcProjRHS(matchfun))

    def matchPointFunction(self, matchfun, r0):
        # NeedUnitTest
        self._solveForRHS(self.calcProjPointfunRHS(matchfun, r0))
        
    def matchErrRMS(self, matchfun):
        """
        Calculate the integral of the RMS error for matchfun and this system
        discretiser DOF set. Note this is approximate because quadrature is
        used for the integration. This is defined as the L2 norm of the error
        divided by the L2 norm of the matching function.
        """
        err_intg_acc = 0.
        match_intg_acc = 0.
        for el in self.disc.elements:
            intg = el.rule.integrateFun
            perm_l, perm_g = el.permutation()
            reconstVals = N.sum(el.physVals()[perm_l]*self.dofArray[perm_g][
                :, N.newaxis, N.newaxis], axis=0)
            evMatchfun = N.array(map(matchfun, el.physEvalPoints()))
            err = reconstVals - evMatchfun
            size = el.size
            # Integrate |reconstructed value - matchfun|^2
            err_intg_acc += intg(N.fromiter((N.dot(x,x) for x in err), err.dtype))*size
            match_intg_acc += intg(
                N.fromiter((N.dot(x,x) for x in evMatchfun), err.dtype))*size
        return N.sqrt(N.abs(err_intg_acc/match_intg_acc))

class PformDiscretiserRHS(PformDiscretiserBase):
    def calcProjRHS(self, matchfun):
        print "Calculating Function Projection RHS"
        RHS = N.zeros_like(self.dofArray)
        for el in self.disc.elements:
            intg = el.rule.integrateFun
            physVals = el.physVals()
            evMatchfun = N.array(map(matchfun, el.physEvalPoints()))
            local_RHS = N.array([intg(N.sum(fn_i*evMatchfun, axis=1))
                                 for fn_i in physVals])*el.size
            elPerm = el.permutation()
            RHS[elPerm[1]] += local_RHS[elPerm[0]]
        return RHS

    def calcProjPointfunRHS(self, matchfun, r0):
        RHS = N.zeros_like(self.dofArray)
        local_RHS, elPerm = self.calcProjPointfunRHS_with_elPerm(matchfun, r0)
        RHS[elPerm[1]] = local_RHS[elPerm[0]]
        return RHS

    def calcProjPointfunRHS_with_elPerm(self, matchfun, r0):
        print "Calculating Point Function Projection RHS"
        funval = matchfun(r0)
        elno, l_coord = PostProc.LocatePoints(self.disc.mesh, N.array([r0]))
        elno = elno[0]
        el = self.disc.elements[elno]
        bf_vals = el.physValsAtPoints(l_coord)[:,0]
        local_RHS = N.array([N.dot(val, funval) for val in bf_vals])
        elPerm = el.permutation()
        return local_RHS, elPerm

class PformDiscretiserDOFs(PformDiscretiserRHS):

    def D(self, target_disc):
        target_disc_dofs = target_disc.newDOFs()
        target_disc_dofs.dofArray[:] = self.matrix.exteriorDerivative(
            target_disc)*self.dofArray
        return target_disc_dofs

    def hodgeStar(self, target_disc):
        # Hodge Star of a n-dimensional p-form results in an n-dimensional
        # (n-p) form. H
        n = 3
        if target_disc.p != n - self.disc.p:
            raise TypeError('Hodge start of a %d-form must result in a %d-form'\
                            % (self.disc.p, n - self.disc.p))
        return self.projOnto(target_disc)

    def projOnto(self, target_disc, extDiffTesting=None):
        """
        Projects this set of DOFs onto the DOFs of target_disc.

        This is equivalent to solving [M_target]{x} = [P]{y} for {x} where

        [M_target] -- target mass matrix
        {x} -- target dof vector
        [P] -- Matrix This discretiser's basis functions tested by the target basis functions
        {y} -- This discretiser's dofs.

        If extDiffTesting is specified as True, this discretiser's basis functions are tested
        with the exterior derivative of target_disc's basisfunctions.
        """
        if not extDiffTesting: P = self.matrix.projectionOnto(target_disc)
        else: P = self.matrix.projectionOnto(target_disc.D())
        M = self.useLU and target_disc.matrix.mass_LU() or target_disc.matrix.mass()
        target_dofs = target_disc.newDOFs()
        target_dofs.dofArray[:] = self.solver.solve_mat_matvec(M, P, self.dofArray)
        return target_dofs
        

class PformDiscretiserDOFs_D(PformDiscretiserRHS): pass
