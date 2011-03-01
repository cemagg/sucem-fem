from __future__ import division

import numpy as N
from scipy import sparse

from NewCode.Utilities import Struct, partial
from NewCode.DiscretisedSystem import BrickCoupledFirstOrderSystemBDirechlet, \
     PMLSystem, BrickCoupledFirstOrderSystemB

class CoupledBond(object):
    def __init__(self, bond_sys, master_sys):
        if hasattr(master_sys, "get_bond_data"):
            bond_data = master_sys.get_bond_data()
            master_E_dofs = bond_data.E_dofs
            master_E_disc = bond_data.E_disc
        else:
            master_E_dofs = master_sys.dofs.E
            master_E_disc = master_sys.discs.E
        bond_B_dofs = bond_sys.dofs.B
        bond_B_disc = bond_sys.discs.B
        self.e_m = master_E_dofs.dofArray
        self.b_b = bond_B_dofs.dofArray
        self.Mb_b = bond_B_disc.matrix.mass()
        self.C_bm = master_E_disc.matrix.partialExteriorDerivative(
            bond_B_disc)


class CoupledBondingSystem(BrickCoupledFirstOrderSystemB):
    def __init__(self, in_self, in_self_vol=None):
        self.bond_surface_fns = {}
        self.in_self = in_self
        if  in_self_vol: self.in_self_vol = in_self_vol
        else: self.in_self_vol = in_self
        
    @staticmethod
    def calc_bonding_freefuns(glob_freefuns, in_self, bond_surface_fns):
        gfreeE = glob_freefuns.E ; gfreeB = glob_freefuns.B
        freeE = lambda ent: gfreeE(ent) and not N.any(
            [bs(ent) for bs in bond_surface_fns]) and in_self(ent)
        freeB = lambda ent: gfreeB(ent) and not N.any(
            [bs(ent) for bs in bond_surface_fns]) and in_self(ent)
        return Struct(E=freeE, B=freeB)
    
    def add_bonded_system(self, system, bond_surface_fn):
        if hasattr(self, 'discs'): raise Exception(
            "Bonded systems can't be added after init_discs() has been called")
        self.bond_surface_fns[system] = bond_surface_fn

    def init_discs(self, *names, **kwargs):
        b_s_fns = self.bond_surface_fns
        kwargs['btype'] = 'cohen98'
        kwargs['BCs'] = self.calc_bonding_freefuns(
            kwargs['BCs'], self.in_self, b_s_fns.values())
        kwargs['volBCs'] = Struct(E=self.in_self_vol, B=self.in_self_vol)
        BrickCoupledFirstOrderSystemB.__init__(self, *names, **kwargs)
        order = self.discOrders['E'][0]
        assert(order == self.discOrders['B'][0])
        self.discs.E.diagonalise()
        self.discs.B.diagonalise()
        bonded_systems = b_s_fns.keys()
        self.bonds = dict((bs, CoupledBond(self, bs))
                          for bs in bonded_systems)

    def get_bond(self, master_sys):
        return self.bonds[master_sys]

    def step_B(self, no_steps=1):
        print 'stepping CoupledBonding B ', no_steps
        dofs = self.dofs
        discs = self.discs
        dt = self.dt
        e = dofs.E.dofArray ; b = dofs.B.dofArray
        try: C = discs.E.matrix.exteriorDerivative(discs.B)
        except TypeError: # It's possible for the bonding system to have 0 E-dofs
            C_zeros = N.zeros_like(b)
            C = Struct(matvec=lambda x: C_zeros)
        bonds = [(bond.C_bm, bond.e_m) for bond in self.bonds.values()]
        for step in xrange(no_steps):
            # d/dt B = -curl(E_self) - curl(E_bonded)
            b -= dt*N.sum(
                [C*e] + [C_bm*e_m for C_bm, e_m in bonds],
                axis=0)
            yield

    def step_E(self, no_steps=1):
        print 'stepping CoupledBonding E ', no_steps
        dofs = self.dofs
        discs = self.discs
        dt = self.dt
        if discs.E.totalDOFs == 0: # bonding system can have 0 E-dofs
            for step in xrange(no_steps):
                self.n += 1
                yield                   # so do nothing
        C = discs.E.matrix.exteriorDerivative(discs.B)
        Me = self.useLU and discs.E.matrix.mass_LU() or discs.E.matrix.mass()
        Me_solve = partial(dofs.E.solver.solve_mat_vec, Me)
        Mb = discs.B.matrix.mass()
        e = dofs.E.dofArray ; b = dofs.B.dofArray
        for step in xrange(no_steps):
            self.n += 1
#             print 'Bonding Step %d/%d, total steps: %d, max E: %f ' % \
#                   (step+1, no_steps, self.n, N.max(N.abs(self.dofs.E.dofArray)))
            # d/dt E = curl(B)
            e += dt*Me_solve(C.T*(Mb*b))
            self.log()
            yield

class BondedSystem(object):
    def set_bond(self, bond):
        self.bond = bond

class CoupledBondedSystem(BondedSystem, BrickCoupledFirstOrderSystemBDirechlet):
    dynamicDirechlet = False            # do direchlet BC DOFs change at each timestep

    def __init__(self, *names, **kwargs):
        kwargs['btype'] = 'cohen98'
        BrickCoupledFirstOrderSystemBDirechlet.__init__(self, *names, **kwargs)
        self.discs.E.diagonalise()
        self.discs.B.diagonalise()

    def set_dyn_direch(self, dyn_direch_fn):
        self.dynamicDirechlet = True
        self.dyn_direch_fn = dyn_direch_fn
        self.no_dynB_allzeros = 0
        self.dynB_started = False
    
    def get_B_driveContribs(self, dt, n):
        if self.dynamicDirechlet: return self.get_B_driveContribs_dyn(dt,n)
        else: return self.get_B_driveContribs_static(dt, n)

    def get_B_driveContribs_static(self, dt, n):
        try: B_drv_dofs, B_drv_weights = self._B_driveContribs_static
        except AttributeError: 
            B_drv_dofs, B_drv_weights = self._B_driveContribs_static = (
                BrickCoupledFirstOrderSystemBDirechlet.get_B_driveContribs(self))
        try: drive_fun = self.drive_fun
        except AttributeError: drive_fun = lambda *x: 0.
        return B_drv_dofs, B_drv_weights*drive_fun(dt, n)

    def get_B_driveContribs_dyn(self, dt, n):
        da = self.direchSys.dofs.E.dofArray
        if self.no_dynB_allzeros <= 1:
            da[:] = self.dyn_direch_fn(dt, n)
        else: print "Skipping dynB calc"
        if N.max(N.abs(da)) == 0:
            if self.dynB_started:
                self.no_dynB_allzeros +=1
        else: self.dynB_started = True
        return BrickCoupledFirstOrderSystemBDirechlet.get_B_driveContribs(self)
    
    def step_B(self, no_steps=1):
        print 'stepping CoupledBonded B ', no_steps
        dofs = self.dofs
        discs = self.discs
        dt = self.dt
        C = discs.E.matrix.exteriorDerivative(discs.B)
        e = dofs.E.dofArray ; b = dofs.B.dofArray
        try: drive_fun = self.drive_fun
        except AttributeError: drive_fun = lambda *x: 0.
        for step in xrange(no_steps):
            B_drv_dofs, B_drv_weights = self.get_B_driveContribs(dt,self.n)
#             print 'BondedStep %d/%d, total steps: %d, max E: %f ' % \
#                   (step+1, no_steps,  self.n, N.max(N.abs(self.dofs.E.dofArray)))
            # d/dt B = -curl(E_self) 
            b -= dt*C*e
            b[B_drv_dofs] -= dt*B_drv_weights
            yield

    def step_E(self, no_steps=1):
        print 'stepping CoupledBonded E ', no_steps
        dofs = self.dofs
        discs = self.discs
        dt = self.dt
        C = discs.E.matrix.exteriorDerivative(discs.B)
        C_bm = self.bond.C_bm ; Mb_b = self.bond.Mb_b
        b_b = self.bond.b_b
        Me = self.useLU and discs.E.matrix.mass_LU() or discs.E.matrix.mass()
        Me_solve = partial(dofs.E.solver.solve_mat_vec, Me)
        Mb = discs.B.matrix.mass()
        E_RHS_drive = self.get_E_RHS_driveContribs()
        E_drv_dofsd, E_drv_weightsd = E_RHS_drive.delta
        E_drv_dofsc, E_drv_weightsc = E_RHS_drive.current        
        try: drive_fun = self.drive_fun
        except AttributeError: drive_fun = lambda *x: 0.
        e = dofs.E.dofArray ; b = dofs.B.dofArray
        for step in xrange(no_steps):
            drv_np1 = drive_fun(dt, self.n+1)
            drv_n = drive_fun(dt, self.n)
            drv_delta = drv_np1 - drv_n
            self.n += 1
            # d/dt E = curl(B_self) + curl(B_bonding)
            RHS = dt*(C.T*(Mb*b) + C_bm.T*(Mb_b*b_b))
            RHS[E_drv_dofsd] -= drv_delta*E_drv_weightsd 
            RHS[E_drv_dofsc] -= dt*drv_n*E_drv_weightsc
            e += Me_solve(RHS)
            self.log()
            yield
    
    def reset_history(self, *names, **kwargs):
        super(CoupledBondedSystem, self).reset_history(*names, **kwargs)
        if self.dynamicDirechlet:
            self.no_dynB_allzeros = 0
            self.dynB_started = 0
            


class PMLBondedSystem(BondedSystem, PMLSystem):
    def step_B(self, no_steps=1):
        print 'stepping PMLBondedSystem', no_steps
        dofs = self.dofs
        discs = self.discs 
        dt = self.dt
        print 'Calculating Curl matrix'
        C = discs.E.matrix.exteriorDerivative(discs.B)
        A_by = self.PMLMats.A_by() ; B_by = self.PMLMats.B_by()
        A_bx = self.PMLMats.A_bx() ; B_bx = self.PMLMats.B_bx()
        A_hz = self.PMLMats.A_hz() ; B_hz = self.PMLMats.B_hz()
        # We are copying the matrix here, should try to access .data directly,
        # but diag matrix works differently to csc/coo/csr
        L_2 = discs.B.matrix.mass().diagonal() # matrix should be diagonal!
        B_drv_dofs, B_drv_weights = self.get_B_driveContribs()
        try: drive_fun = self.drive_fun
        except AttributeError: drive_fun = lambda *x: 0.
        b = dofs.B.dofArray ; h = dofs.h.dofArray
        e = dofs.E.dofArray ; d = dofs.d.dofArray
        for step in xrange(no_steps):
            drv_n = drive_fun(dt, self.n)
            # Update fake B using curl of E
            next_b = 1/A_by*(B_by*b - L_2*(C*e))
            next_b[B_drv_dofs] -= drv_n/A_by[B_drv_dofs]*L_2[B_drv_dofs] \
                                  *B_drv_weights
            # Update dual-mesh H using fake B
            h[:] = 1/A_hz*(A_bx*next_b - B_bx*b + B_hz*h)
            b[:] = next_b    
            del(next_b)
            yield
            
    def step_E(self, no_steps=1):
        E_RHS_drive = self.get_E_RHS_driveContribs()
        E_drv_dofs, E_drv_weights = E_RHS_drive.current
        try: drive_fun = self.drive_fun
        except AttributeError: drive_fun = lambda *x: 0.
        C = self.discs.E.matrix.exteriorDerivative(self.discs.B)
        C_bm = self.bond.C_bm ; Mb_b = self.bond.Mb_b
        b_b = self.bond.b_b
        A_dy = self.PMLMats.A_dy() ; B_dy = self.PMLMats.B_dy()
        A_dx = self.PMLMats.A_dx() ; B_dx = self.PMLMats.B_dx()
        A_ez = self.PMLMats.A_ez() ; B_ez = self.PMLMats.B_ez()
        h = self.dofs.h.dofArray
        e = self.dofs.E.dofArray ; d = self.dofs.d.dofArray
        dt = self.dt
        for step in xrange(no_steps):
            drv_n = drive_fun(dt, self.n)
            self.n += 1
#             print 'PML Step %d/%d, drv_fun: %f, total steps: %d, max E: %f ' % \
#                   (step+1, no_steps, drv_n, self.n, N.max(N.abs(self.dofs.E.dofArray)))
            # Update fake d using curl of H and -J
            next_d = 1/A_dy*(C.T*h + B_dy*d + C_bm.T*(Mb_b*b_b))
            next_d[E_drv_dofs] -= drv_n/A_dy[E_drv_dofs]*E_drv_weights
            # Update E using fake D
            e[:] = 1/A_ez*(A_dx*next_d - B_dx*d + B_ez*e)
            d[:] = next_d    
            del(next_d)
            self.log()            
            yield
            
