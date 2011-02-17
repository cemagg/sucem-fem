from __future__ import division

import numpy as N

class E_intg(object):
    """
    Calculate the integral (n x w_i)exp(j*k*unit(r).r') over Gamma_C where 
    w_i is the i'th E-field basis fun,
    Gamma_C the NTFF collection surface,
    k the wavenumber
    unit(r) the radial unit vector to the observation point
    r' the local source coordinate on Gamma_C.

    Gamma_C and the w_i's are specified by E_subdisc
    """

    def __init__(self, E_subdisc):
        self.subdisc = E_subdisc

    def calc(self, r_hat, k, fft_dofs):
        assert(N.allclose(N.linalg.norm(r_hat), 1.))
        spdim = 3
        acc = N.zeros(spdim, dtype=N.complex128)
        for el in self.subdisc.elements:
            r_prime = evpts = el.physEvalPoints()
            ejkrr = N.exp(1j*k*N.sum(r_hat*r_prime, axis=1))
            intg = el.rule.integrateFun
            perm_l, perm_g = el.permutation()
            dof_vals = fft_dofs[perm_g]
            pvals = el.physVals()
            fnvals = pvals*ejkrr[N.newaxis,:,N.newaxis]
            intg_vals = [intg(fnv) for fnv in  fnvals]
            acc += N.sum(intg_vals*dof_vals[:,N.newaxis], axis=0)*el.area()

        return acc
    
