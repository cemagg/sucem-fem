from __future__ import division

import numpy as np

def S11(Zl, Z0):
    """Calculate reflection coeficient of load

    @param Zl: Load impedance
    @param Z0: Characteristic impedance
    """
    Zl = np.asarray(Zl)
    return (Zl - Z0)/(Zl + Z0)
