from __future__ import division

import numpy as N
from scipy import sparse
from scipy.sparse.linalg import dsolve as linsolve
from scipy.io import numpyio
from NewCode.Utilities import CacheLast, Struct, max_or_scalar
from NewCode import PostProc, Integration
from NewCode.MatrixUtils import Matrices
from NewCode.DifferentialForm import Discretiser, BrickDiscretiser, allconstrained
from NewCode.DifferentialForm import constrained_on_boundary, DiscretiserEntities

class EigenSystem(object):
    """
    Setup eigensystem for a p-form discretisation. Subclasses need to define p.
    """
    p = None                            # Use p-form to discretise system
    DiscretiserModule = Discretiser

    def __init__(self, mesh, order, mixed=True, BC=constrained_on_boundary,
                 volBC=None, btype=None,  **kwargs):
        self.order = order
        self.mesh = mesh

        if not btype:
            self.disc = self.DiscretiserModule.setup_PformDiscretiser(
                mesh, self.p, order, mixed, BC, volBC, **kwargs)
        else:
            self.disc = self.DiscretiserModule.setup_PformDiscretiser(
                mesh, self.p, order, mixed, BC, volBC, btype=btype,
                **kwargs)
        self.btype = btype
        self.massMatrix = self.disc.matrix.mass
        self.stiffnessMatrix = self.disc.matrix.stiffness
    

class VectorWaveEigen(EigenSystem):
    p = 1                               # E/H fields are 1-forms

class AudioEigen(EigenSystem):
    p = 2                               # Audio velocity is 2-form

class TimeDomainSystem(object):
    def zero(self):
        """ Set all DOFs to 0 """
        for dof_name,dof in self.dofs.iteritems(): dof.dofArray[:] = 0

    def reset_history(self, zero_dofs=True):
        if zero_dofs: self.zero()
        self.n = 0
        self.initLogs()

class CurlCurlNewmark(TimeDomainSystem,EigenSystem):
    p = 1
    useLU = False
    faceIntegrationOrder = 8
    faceIntegrator = Integration.TriIntegrator
    def __init__(self, mesh, beta=0.25001, useQ=False, *names, **kwargs):
        super(CurlCurlNewmark, self).__init__(mesh, *names, **kwargs)
        self.beta = beta
        self.dofs = self.disc.newDOFs()
        self.dofArrays = {-1:self.dofs.dofArray.copy(),
                          0:self.dofs.dofArray.copy(),
                          1:self.dofs.dofArray}
        self.solver = self.dofs.solver
        self.hasQ = useQ
        if useQ:
            self.boundarySurfaceMatrix = self.disc.matrix.boundarySurfaceMatrix
            self.disc.elements[:].setFaceIntegrationRule(
                self.faceIntegrator(self.faceIntegrationOrder))
        self.reset_history(zero_dofs=False)
        
    def setTimestep(self,dt):
        self.dt = dt
        self.solveMat = 1/(dt**2)*self.massMatrix() + self.beta*self.stiffnessMatrix()
        if self.hasQ:
            self.solveMat = self.solveMat + 1/2/dt*self.boundarySurfaceMatrix()
        if self.useLU:
            self.solveMat_orig = self.solveMat
            self.solveMat = Struct(solve=linsolve.factorized(self.solveMat))
        
    def setDriveDOFs(self, dofnos, weights, drv_fun):
        assert(N.isscalar(weights) or len(weights) == len(dofnos))
        assert(max_or_scalar(dofnos) < self.disc.totalDOFs)
        self.drive = (dofnos, weights, drv_fun)

    set_driveDOFs = setDriveDOFs

    def initLogs(self):
        self.loggedDOFs = {}
        self.loggedReconstructed = []

    def addLogger(self, dofnos, divisor=1):
        """
        Log values of DOF # dofnos values for discretiser discName every divisor timesteps
        """
        assert max(dofnos) < self.disc.totalDOFs
        self.loggedDOFs[tuple(dofnos)] = Struct(div=divisor, vals=[])

    def addReconstructedLogger(self, log_elnos, log_el_coords, divisor=1):
        self.loggedReconstructed.append(Struct(
            reconparms=list(self.dofs.calc_reconparms(
            log_elnos, log_el_coords)), div=divisor, vals=[]))

    def log(self):
        for dofNos, log in self.loggedDOFs.items():
            if self.n % log.div == 0:
                log.vals.append(self.dofs.dofArray[list(dofNos)])
        for rLog in self.loggedReconstructed:
            if self.n % rLog.div == 0:
                dofs = self.dofs
                rLog.vals.append(dofs.recon_fromparms(
                    rLog.reconparms, dofs.dofArray))

    def step(self, no_steps=1):
        M = self.massMatrix()
        S = self.stiffnessMatrix()
        if self.hasQ: Q = self.boundarySurfaceMatrix()
        else: Q = Matrices.sparseType(M.shape, nzmax=0)
        dt = self.dt
        b = self.beta
        drv_dofnos, drv_weights, drv_fun = self.drive
        for step in xrange(no_steps):
            self.n += 1
            self.dofArrays[-1][:] = self.dofArrays[0]
            self.dofArrays[0][:] = self.dofArrays[1]
            dm1, d0 = self.dofArrays[-1], self.dofArrays[0]
            RHS = (2/(dt**2)*M.matvec(d0) - (1-2*b)*S.matvec(d0)) \
                  - (1/(dt**2)*M.matvec(dm1) - 1/2/dt*Q.matvec(dm1)
                     + b*S.matvec(dm1))
            drvs = [drv_fun(dt, self.n - i) for i in range(3)]
            RHS[drv_dofnos] -= drv_weights*(
                b*drvs[0] + (1-2*b)*drvs[1] + b*drvs[2])
            print 'Step %d/%d, drv_fun: %f, total steps: %d, max: %f' % \
                    (step+1, no_steps, drvs[0], self.n, N.max(N.abs(self.dofArrays[1])))
            self.dofArrays[1][:] = self.solver.solve_mat_vec(
                self.solveMat, RHS)
            self.log()

    def zero(self):
        for da in self.dofArrays.values(): da[:]=0

class CurlCurlNewmarkDirechlet(CurlCurlNewmark):
    def setDirechBCs(self, DirechBCs, DirechVolBCs=allconstrained):
        self.direchSys = self.__class__(
            self.mesh, order=self.order, BC=DirechBCs, volBC=DirechVolBCs,
            useQ=self.hasQ, btype=self.btype)
    def setTimestep(self,dt):
        super(CurlCurlNewmarkDirechlet, self).setTimestep(dt)
        disc_f, disc_p = self.disc, self.direchSys.disc
        Mp = disc_p.matrix.projectionOnto(disc_f)
        Sp = disc_p.D().matrix.projectionOnto(disc_f.D())
        self.direchletSolveMat = 1/(dt**2)*Mp + self.beta*Sp

    def step(self, no_steps=1):
        Mf = self.massMatrix()
        Sf = self.stiffnessMatrix()
        if self.hasQ: raise NotImplementedError
        Qf = Matrices.sparseType(Mf.shape, nzmax=0)
        disc_f, disc_p = self.disc, self.direchSys.disc
        Mp = disc_p.matrix.projectionOnto(disc_f)
        Sp = disc_p.D().matrix.projectionOnto(disc_f.D())
        Qp = Matrices.sparseType(Mp.shape, nzmax=0)
        direch_dofs = self.direchSys.dofs.dofArray
        direch_solveMat = self.direchletSolveMat
        dt = self.dt
        b = self.beta
        drv_fun = self.drive_fun
        for step in xrange(no_steps):
            self.n += 1
            self.dofArrays[-1][:] = self.dofArrays[0]
            self.dofArrays[0][:] = self.dofArrays[1]
            dm1, d0 = self.dofArrays[-1], self.dofArrays[0]
            RHS_f = (2/(dt**2)*Mf.matvec(d0) - (1-2*b)*Sf.matvec(d0)) \
                    - (1/(dt**2)*Mf.matvec(dm1) - 1/2/dt*Qf.matvec(dm1)
                       + b*Sf.matvec(dm1))
            drvs = [drv_fun(dt, self.n - i) for i in range(3)]
            # drvs[0,1,2] -> drv(n+1), drv(n), drv(n-1)
            dp1_p, d0_p, dm1_p = [drv*direch_dofs for drv in drvs]
            RHS_d = (2/(dt**2)*Mp.matvec(d0_p) - (1-2*b)*Sp.matvec(d0_p)) \
                    - (1/(dt**2)*Mp.matvec(dm1_p) - 1/2/dt*Qp.matvec(dm1_p)
                       + b*Sp.matvec(dm1_p)) \
                    - direch_solveMat.matvec(dp1_p)
            RHS = RHS_f + RHS_d
            print 'Newmark step %d/%d, drv_fun: %f, total steps: %d, max: %f' % \
                    (step+1, no_steps, drvs[0], self.n, N.max(N.abs(self.dofArrays[1])))
            self.dofArrays[1][:] = self.solver.solve_mat_vec(
                self.solveMat, RHS)
            self.log()

class CurlCurlNewmarkDirechletCoupled(CurlCurlNewmarkDirechlet):

    def setDirechSystem(self, direch_sys):
        self.direchSys = direch_sys

    def step(self, no_steps=1):
        Mf = self.massMatrix()
        Sf = self.stiffnessMatrix()
        if self.hasQ: raise NotImplementedError
        Qf = Matrices.sparseType(Mf.shape, nzmax=0)
        disc_f, disc_p = self.disc, self.direchSys.disc
        Mp = disc_p.matrix.projectionOnto(disc_f)
        Sp = disc_p.D().matrix.projectionOnto(disc_f.D())
        Qp = Matrices.sparseType(Mp.shape, nzmax=0)
        direch_solveMat = self.direchletSolveMat
        dt = self.dt
        b = self.beta
        for step in xrange(no_steps):
            self.n += 1
            self.dofArrays[-1][:] = self.dofArrays[0]
            self.dofArrays[0][:] = self.dofArrays[1]
            dm1, d0 = self.dofArrays[-1], self.dofArrays[0]
            RHS_f = (2/(dt**2)*Mf.matvec(d0) - (1-2*b)*Sf.matvec(d0)) \
                    - (1/(dt**2)*Mf.matvec(dm1) - 1/2/dt*Qf.matvec(dm1)
                       + b*Sf.matvec(dm1))
            # Direchlet DOF values at t = n+1, n, n-1
            dp1_p, d0_p, dm1_p = self.direchSys.dof_mem
            RHS_d = (2/(dt**2)*Mp.matvec(d0_p) - (1-2*b)*Sp.matvec(d0_p)) \
                    - (1/(dt**2)*Mp.matvec(dm1_p) - 1/2/dt*Qp.matvec(dm1_p)
                       + b*Sp.matvec(dm1_p)) \
                    - direch_solveMat.matvec(dp1_p)
            RHS = RHS_f + RHS_d
            print 'Newmark step %d/%d, total steps: %d, max_RHS_d: %f, max_RHS_f: %f, max: %f' % \
                    (step+1, no_steps, self.n, N.max(N.abs(RHS_d)), N.max(N.abs(RHS_f)),
                     N.max(N.abs(self.dofArrays[1])))
            self.dofArrays[1][:] = self.solver.solve_mat_vec(
                self.solveMat, RHS)
            self.log()
 

class CoupledFirstOrderSystem(TimeDomainSystem):
    useLU = False
    DiscretiserModule = Discretiser
    discForms = {'E':1, 'H':1, 'D':2, 'B':2} # 'disc_name': disc_form
    discNames = tuple(discForms.keys())
    def __init__(self, mesh, BCs=None, volBCs=None, disc_orders=None, init_disc=True,
                 btype=None, el_in_disc=None):
        self.discs = Struct()
        self.dofs = Struct()
        self.initLogs()
        self.mesh = mesh
        if disc_orders == None:
            disc_orders = {'E':(1,True), 'H':(1,True), 'D':(1,True), 'B':(1,True)}
        if volBCs == None:
            volBCs = dict(E=None, H=None, D=None, B=None)
        self.discOrders = disc_orders
        self.btype = btype
        if init_disc:
            for disc_name in self.discForms:
                self.initDiscretiser(
                    disc_name, BCs[disc_name], volBCs[disc_name], btype=btype,
                    el_in_disc=el_in_disc)
        self.reset_history(zero_dofs=False)

    def setTimestep(self,dt):
        for dof in self.dofs.values(): dof.useLU = self.useLU
        self.dt = dt

    def initDiscretiser(self, disc_name, BC, volBC=None, btype=None, el_in_disc=None):
        order, mixed = self.discOrders[disc_name]
        p = self.discForms[disc_name]
        mesh = self.mesh
        self.btype = btype
        if not btype:
            self.discs[disc_name] = self.DiscretiserModule.setup_PformDiscretiser(
                mesh, p, order, mixed, BC, volBC, el_in_disc=el_in_disc)
        else:
            self.discs[disc_name] = self.DiscretiserModule.setup_PformDiscretiser(
                mesh, p, order, mixed, BC, volBC, btype=btype, el_in_disc=el_in_disc)
        self.dofs[disc_name] = self.discs[disc_name].newDOFs()

    def setDriveDOFs_J(self, dofnos, weights, drv_fun):
        assert(N.isscalar(weights) or len(weights) == len(dofnos))
        assert(max_or_scalar(dofnos) < self.discs.D.totalDOFs)
        self.drive_J = (dofnos, weights, drv_fun)

    setDriveDOFs = setDriveDOFs_J

    def step(self, no_steps=1):
        print 'stepping ', no_steps
        dofs = self.dofs
        dt = self.dt
        for step in xrange(no_steps):
            dofnos, weights, drv_fun = self.drive_J
            drv = drv_fun(dt, self.n)
            print 'Step %d/%d, drv_fun: %f, total steps: %d, max E: %f ' % \
                  (step+1, no_steps, drv, self.n, N.max(N.abs(self.dofs.E.dofArray)))
            dofs.D.dofArray[:] += dt*dofs.H.D(dofs.D.disc).dofArray
            dofs.D.dofArray[dofnos] -= dt*weights*drv
            dofs.E = dofs.D.hodgeStar(dofs.E.disc)
            self.log()
            dofs.B.dofArray[:] -= dt*dofs.E.D(dofs.B.disc).dofArray
            dofs.H = dofs.B.hodgeStar(dofs.H.disc)
            self.n += 1
            
    def log(self):
        for discName, logs in self.loggedDOFs.items():
            for dofNos, log in logs.items():
                if self.n % log.div == 0:
                    log.vals.append(self.dofs[discName].dofArray[list(dofNos)])
            for fLog in self.loggedFileDOFs[discName]:
                if self.n % fLog.div == 0:
                    numpyio.fwrite(fLog.file, len(fLog.dofnos),
                                   self.dofs[discName].dofArray[fLog.dofnos])
                    fLog.len += 1
            for rLog in self.loggedReconstructed[discName]:
                if self.n % rLog.div == 0:
                    dofs = self.dofs[discName]
                    rLog.vals.append(dofs.recon_fromparms(
                        rLog.reconparms, dofs.dofArray))
                    
    def initLogs(self):
        self.loggedDOFs = Struct((k,{}) for k in self.discNames)
        self.loggedReconstructed = Struct((k,[]) for k in self.discNames)
        self.loggedFileDOFs = Struct((k,[]) for k in self.discNames)

    def addLogger(self, discName, dofnos, divisor=1):
        """
        Log values of DOF # dofnos values for discretiser discName every divisor timesteps
        """

        assert max(dofnos) < self.discs[discName].totalDOFs
        self.loggedDOFs[discName][tuple(dofnos)] = Struct(div=divisor, vals=[])

    def addFileLogger(self, discName, dofnos, fileobj, divisor=1):
        """
        Log values of DOF # dofnos values for discretiser discName every divisor timesteps
        """

        assert max(dofnos) < self.discs[discName].totalDOFs
        self.loggedFileDOFs[discName].append(Struct(
            file=fileobj, div=divisor, dofnos=dofnos, len=0))
        type_char = self.dofs[discName].dofArray.dtype.char
        fileobj.write(type_char + ' ' + str(len(dofnos)) + ' \n')

    def addReconstructedLogger(self, discName, log_elnos, log_el_coords, divisor=1):
        self.loggedReconstructed[discName].append(Struct(
            reconparms=list(self.dofs[discName].calc_reconparms(
            log_elnos, log_el_coords)), div=divisor, vals=[]))


class CoupledFirstOrderSystemHardSource(CoupledFirstOrderSystem):

    def setSourceDOFs(self, weights, drv_fun):
        self.sourceDOFs = (weights, drv_fun)

    def setDirechBCs(self, DirechBCs):
        volBCs=dict(E=allconstrained, B=allconstrained,
                    H=allconstrained, D=allconstrained)
        self.direchSys = CoupledFirstOrderSystem(
            self.mesh, DirechBCs, volBCs, disc_orders=self.discOrders,
            btype=self.btype)

    def step(self, no_steps=1):
        print 'stepping ', no_steps
        dofs = self.dofs
        direch_dofs = self.direchSys.dofs
        dt = self.dt
        weights, drv_fun = self.sourceDOFs
        for step in xrange(no_steps):
            self.n += 1
            drv = drv_fun(dt, self.n)
            print 'Step %d/%d, drv_fun: %f, total steps: %d, max E: %f ' % \
                  (step+1, no_steps, drv, self.n, N.max(N.abs(self.dofs.E.dofArray)))
            dofs.D.dofArray[:] += dt*dofs.H.D(dofs.D.disc).dofArray
            dofs.E = dofs.D.hodgeStar(dofs.E.disc)
            dofs.B.dofArray[:] -= dt*dofs.E.D(dofs.B.disc).dofArray
            direch_dofs.E.dofArray[:] = weights*drv
            dofs.B.dofArray[:] -= dt*direch_dofs.E.D(dofs.B.disc).dofArray
            dofs.H = dofs.B.hodgeStar(dofs.H.disc)
            self.log()
    

class CoupledFirstOrderSystemB(CoupledFirstOrderSystem):
    discForms = {'E':1, 'B':2}
    discNames = tuple(discForms.keys())

    def setDriveDOFs_J(self, dofnos, weights, drv_fun):
        assert(N.isscalar(weights) or len(weights) == len(dofnos))
        assert(max_or_scalar(dofnos) < self.discs.E.totalDOFs)
        self._drive_J = (dofnos, weights, drv_fun)
        self.drive_fun = drv_fun

    setDriveDOFs = setDriveDOFs_J
        
    def drive_J(self):
        try: dofnos, weights, drv_fun = self._drive_J
        except AttributeError:
            dofnos, weights, drv_fun = [0], [0], lambda *x: 0
        return dofnos, weights, drv_fun

    def step(self, no_steps=1):
        print 'stepping ', no_steps
        dofs = self.dofs
        discs = self.discs
        dt = self.dt
        Me = discs.E.matrix.mass_LU() if self.useLU else discs.E.matrix.mass()
        Cbe = discs.B.matrix.projectionOnto(discs.E.D())
        dofnos, weights, drv_fun = self.drive_J()
        for step in xrange(no_steps):
            self.n += 1
            drv = drv_fun(dt, self.n)
            # d/dt B = -curl(E)
            dofs.B.dofArray[:] -= dt*dofs.E.D(dofs.B.disc).dofArray
            # d/dt E = curl(B) - J

            RHS = dt*Cbe.matvec(dofs.B.dofArray)
            RHS[dofnos] -= dt*weights*drv
            print 'Step %d/%d, drv_fun: %f, total steps: %d, max E: %f ' % \
                  (step+1, no_steps, drv, self.n, N.max(N.abs(self.dofs.E.dofArray)))
            dofs.E.dofArray[:] += dofs.E.solver.solve_mat_vec(Me, RHS)
            self.log()

class CoupledFirstOrderSystemBDirechlet(CoupledFirstOrderSystemB):
        
    def setDirechBCs(
        self, DirechBCs, volBCs=dict(E=allconstrained, B=allconstrained),
        el_in_disc=None):
        self.direchSys = self.__class__(
            self.mesh, DirechBCs, volBCs, self.discOrders, btype=self.btype,
            el_in_disc=el_in_disc)

    def get_E_RHS_driveContribs(self):
        contribs = Struct(current=([],N.array([])), delta=([],N.array([])))
        try:        
            dofnos_J, weights_J, drv_fun = self._drive_J
            contribs.current=(dofnos_J, weights_J)
        except AttributeError: pass
        if hasattr(self, 'direchSys'):
            dofs_p = self.direchSys.dofs
            discs_f, discs_p = self.discs, self.direchSys.discs        
            Me_p = discs_p.E.matrix.projectionOnto(discs_f.E)
            cont = Me_p.matvec(dofs_p.E.dofArray)
            dofnos_direch = N.arange(
                Me_p.shape[0], dtype=N.int32)[N.abs(cont) > 0].copy()
            weights_direch = cont[dofnos_direch]
            contribs.delta = (dofnos_direch, weights_direch)

        return contribs

    def get_B_driveContribs(self):
        contribs = (N.array([], N.int32), N.array([]))
        if hasattr(self, 'direchSys'):
            dofs_p = self.direchSys.dofs
            discs_f, discs_p = self.discs, self.direchSys.discs        
            Ceb_p = discs_p.E.matrix.exteriorDerivative(discs_f.B)
            cont = Ceb_p.matvec(dofs_p.E.dofArray)
            dofnos_direch = N.arange(
                Ceb_p.shape[0], dtype=N.int32)[N.abs(cont) > 0].copy()
            weights_direch = cont[dofnos_direch]
            contribs = (dofnos_direch, weights_direch)
        return contribs
                                   
    def step(self, no_steps=1):
        print 'stepping coupled', no_steps
        dofs = self.dofs
        discs_f = self.discs
        dt = self.dt
        print 'Getting Curl matrix'
        Ceb_f = discs_f.E.matrix.exteriorDerivative(discs_f.B)
        print 'done'
        Me_f = self.useLU and discs_f.E.matrix.mass_LU() or discs_f.E.matrix.mass()
        #Cbe = discs_f.B.matrix.projectionOnto(discs_f.E.D())
        Mb = discs_f.B.matrix.mass()
        E_RHS_drive = self.get_E_RHS_driveContribs()
        E_drv_dofsd, E_drv_weightsd = E_RHS_drive.delta
        E_drv_dofsc, E_drv_weightsc = E_RHS_drive.current        
        B_drv_dofs, B_drv_weights = self.get_B_driveContribs()
        
        for step in xrange(no_steps):
            drv_np1 = self.drive_fun(dt, self.n+1)
            drv_n = self.drive_fun(dt, self.n)
            drv_delta = drv_np1 - drv_n
            self.n += 1        
            print 'Step %d/%d, drv_fun: %f, total steps: %d, max E: %f ' % \
                  (step+1, no_steps, drv_np1, self.n, N.max(N.abs(self.dofs.E.dofArray)))
            # d/dt B = -curl(E)
            dofs.B.dofArray[:] -= dt*(Ceb_f.matvec(dofs.E.dofArray))
            dofs.B.dofArray[B_drv_dofs] -= dt*drv_n*B_drv_weights
            # d/dt E = curl(B) - J
            RHS = (dt*Ceb_f.T.matvec(Mb.matvec(dofs.B.dofArray)))
            RHS[E_drv_dofsd] -= drv_delta*E_drv_weightsd 
            RHS[E_drv_dofsc] -= dt*drv_n*E_drv_weightsc
            dofs.E.dofArray[:] += dofs.E.solver.solve_mat_vec(Me_f, RHS)
            self.log()
            
class BrickCoupledFirstOrderSystemB(CoupledFirstOrderSystemB):
    DiscretiserModule = BrickDiscretiser
    
class BrickCoupledFirstOrderSystemBDirechlet(CoupledFirstOrderSystemBDirechlet):
    DiscretiserModule = BrickDiscretiser    

class BrickCurlCurlNewmark(CurlCurlNewmark):
    DiscretiserModule = BrickDiscretiser
    faceIntegrationOrder = 4
    faceIntegrator = Integration.QuadIntegrator

class BrickCurlCurlNewmarkDirechlet(CurlCurlNewmarkDirechlet):
    DiscretiserModule = BrickDiscretiser

class BrickVectorWaveEigen(VectorWaveEigen):
    DiscretiserModule = BrickDiscretiser

class BrickAudioEigen(AudioEigen):
    DiscretiserModule = BrickDiscretiser

from NewCode.DifferentialForm.PMLMatrices import PMLMatrices

class PMLSystem(BrickCoupledFirstOrderSystemBDirechlet):
    def __init__(self, *names, **kwargs):
        kwargs['btype'] = 'cohen98'
        super(PMLSystem, self).__init__(*names, **kwargs)
        self.dofs.d = Struct(dofArray=N.zeros_like(self.dofs.E.dofArray))
        self.dofs.h = Struct(dofArray=N.zeros_like(self.dofs.B.dofArray))
        order = self.discOrders['E'][0]
        assert(order == self.discOrders['B'][0])
        self.discs.E.diagonalise()
        self.discs.B.diagonalise()
        
    def set_sigmas(self, sigma_fns):
        self.PMLMats = PMLMatrices(self.discs, sigma_fns)

    def setTimestep(self,dt):
        super(PMLSystem, self).setTimestep(dt)
        self.PMLMats.set_dt(dt)
        
    def get_E_RHS_driveContribs(self):
        contribs = Struct(current=([],N.array([])), delta=([],N.array([])))
        try:        
            dofnos_J, weights_J, drv_fun = self._drive_J
            contribs.current=(dofnos_J, weights_J)
        except AttributeError: pass
        return contribs
    
    def step(self, no_steps=1):
        print 'stepping ', no_steps
        dofs = self.dofs
        discs = self.discs 
        dt = self.dt
        print 'Calculating Curl matrix'
        C = discs.E.matrix.exteriorDerivative(discs.B)
        A_dy = self.PMLMats.A_dy() ; B_dy = self.PMLMats.B_dy()
        A_dx = self.PMLMats.A_dx() ; B_dx = self.PMLMats.B_dx()
        A_ez = self.PMLMats.A_ez() ; B_ez = self.PMLMats.B_ez()
        A_by = self.PMLMats.A_by() ; B_by = self.PMLMats.B_by()
        A_bx = self.PMLMats.A_bx() ; B_bx = self.PMLMats.B_bx()
        A_hz = self.PMLMats.A_hz() ; B_hz = self.PMLMats.B_hz()
        # We are copying the matrix here, should try to access .data directly,
        # but diag matrix works differently to csc/coo/csr
        L_2 = discs.B.matrix.mass().diagonal() # matrix should be diagonal!
        E_RHS_drive = self.get_E_RHS_driveContribs()
        E_drv_dofs, E_drv_weights = E_RHS_drive.current
        B_drv_dofs, B_drv_weights = self.get_B_driveContribs()
        
        b = dofs.B.dofArray ; h = dofs.h.dofArray
        e = dofs.E.dofArray ; d = dofs.d.dofArray
        for step in xrange(no_steps):
            drv_n = self.drive_fun(dt, self.n)
            self.n += 1
            print 'Step %d/%d, drv_fun: %f, total steps: %d, max E: %f ' % \
                  (step+1, no_steps, drv_n, self.n, N.max(N.abs(self.dofs.E.dofArray)))
            # Update fake B using curl of E
            next_b = 1/A_by*(B_by*b - L_2*C.matvec(e))
            next_b[B_drv_dofs] -= drv_n/A_by[B_drv_dofs]*L_2[B_drv_dofs] \
                                  *B_drv_weights
            # Update dual-mesh H using fake B
            h[:] = 1/A_hz*(A_bx*next_b - B_bx*b + B_hz*h)
            b[:] = next_b    
            del(next_b)
            # Update fake d using curl of H and -J
            next_d = 1/A_dy*(C.T.matvec(h) + B_dy*d)
            next_d[E_drv_dofs] -= drv_n/A_dy[E_drv_dofs]*E_drv_weights
            # Update E using fake D
            e[:] = 1/A_ez*(A_dx*next_d - B_dx*d + B_ez*e)
            d[:] = next_d    
            del(next_d)
            self.log()            
