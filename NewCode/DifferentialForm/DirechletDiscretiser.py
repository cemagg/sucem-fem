from __future__ import division

import numpy as N
from NewCode.DifferentialForm import Discretiser
from NewCode.Exceptions import AllZero

class DirechletDiscretiser(Discretiser.PformDiscretiser):
    def __init__(self, freeDisc, geomEntities, Permuter, *names, **kwargs):
        fd = freeDisc
        super(DirechletDiscretiser, self).__init__(
            fd.p, fd.mesh, geomEntities, Permuter, fd.basisSet,
            *names, **kwargs)
        
