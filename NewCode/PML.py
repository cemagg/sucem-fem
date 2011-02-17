from __future__ import division

import numpy as N

def calc_PML_sigma_max(h, PML_m=3, nabla=1, sigma_factor=1):
    return 0.8*(PML_m+1)/h/nabla*sigma_factor

class PML_sigma_fn_1DR(object):
    def __init__(self, sigma_max, r_val, l_PML, m=3):
        self.sigma_max = sigma_max
        self.r_val = r_val
        self.m = m
        self.l_PML = l_PML

    def __call__(self, x):
        return N.where(
            x >= self.r_val,
            self.sigma_max*((x-self.r_val)/self.l_PML)**self.m,
            0.)
    
    
class PML_sigma_fn_1DL(object):
    def __init__(self, sigma_max, l_val, l_PML, m=3):
        self.sigma_max = sigma_max
        self.l_val = l_val
        self.m = m
        self.l_PML = l_PML

    def __call__(self, x):
        return N.where(
            x <= self.l_val,
            self.sigma_max*((self.l_val-x)/self.l_PML)**self.m,
            0.)

class PML_sigma_fn_1DLR(object):
    def __init__(self, sigma_max, l_val, r_val, l_PML, m=3):
        self.sigma_max = sigma_max
        self.l_val = l_val
        self.r_val = r_val
        assert (l_val < r_val)
        self.m = m
        self.l_PML = l_PML
        self.l_fn = PML_sigma_fn_1DL(sigma_max, l_val, l_PML, m)
        self.r_fn = PML_sigma_fn_1DR(sigma_max, r_val, l_PML, m)

    def __call__(self, x):
        return self.l_fn(x) + self.r_fn(x)

def fn_1D2fn_3D(fn_1D, coordno):
    def _fn_3D(r):
        return fn_1D(r.T[coordno])
    return _fn_3D


