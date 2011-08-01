from __future__ import division

import numpy as np
import pickle
import pylab
import sys
mpath = '../../'
sys.path.append(mpath)
from sucemfem.PostProcessing.circuit import S11
#sys.path.remove(mpath)

fname = 'patch_o-2_gmres_Z.pickle'
res = pickle.load(open(fname))

Z0 = 50.
Z = res['Zs']
freqs = res['freqs']
S11 = S11(Z, Z0)

pylab.figure(1)
pylab.hold(0)
pylab.plot(freqs, np.real(Z), label='Re')
pylab.hold(1)
pylab.plot(freqs, np.imag(Z), label='Im')
pylab.grid(1)
pylab.legend(loc=0)

pylab.figure(2)
pylab.hold(0)
pylab.plot(freqs, 20*np.log10(np.abs(S11)))
pylab.hold(1)
pylab.grid(1)

