from __future__ import division
"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from itertools import izip, chain
import os, sys
import pickle
import random
import numpy as N
import scipy
from scipy import sparse, linalg
from numpy.testing import assert_almost_equal
#
# Local Imports
#
import NewCode
from NewCode import Utilities, ProxyList
from NewCode.Utilities import Struct,  close_to_point
from NewCode.GeomGen.Hybrid import SurroundedRect

tet_geom_minsize = N.array([0.25,0.25,0.25])
freespace_minsize = N.array([0.5, 0.5, 0.5])

# tet_geom_minsize = N.array([0.5,0.5,0.5])
# freespace_minsize = N.array([0.75, 0.75, 0.75])

#tet_geom_minsize = N.array([2,2,2.])
#freespace_minsize = N.array([2.75, 2.75, 2.75])


h = 1/4.

geom = SurroundedRect(tet_geom_minsize, freespace_minsize, h)
geom.init_background_mesh()
geom.init_dead_background_element_set()
geom.init_geom()

tst_geof = file('+tst_dipole.geo', 'w')
geom.write_geom(tst_geof)
tst_geof.close()
