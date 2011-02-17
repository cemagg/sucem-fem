from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal
from scipy.sparse.linalg import dsolve as linsolve

from NewCode.tests.TestMeshes import TwoTets
from NewCode import Mesh
from NewCode import SystemMatrix
from NewCode.MatrixUtils import Matrices
from NewCode.DifferentialForm import Discretiser
from NewCode.DifferentialForm import DiscretiserMatrices

class test_DiscretiserMatrices(NumpyTestCase):
    def setUp(self):
        self.mesh = Mesh.Mesh(TwoTets.listmesh)
        self.disc1, self.disc2 = map(
            Discretiser.setup_PformDiscretiser, (self.mesh,)*2, (1,2))

    def _test_matrix(self, inst, desired, attr, *discs):
        testmat = getattr(inst, attr)(*discs)
        assert_equal(testmat.todense(), desired)
        self.assert_(testmat is getattr(inst,attr)(*discs))
        getattr(inst, attr).clearCache()
        testmat2 = getattr(inst, attr)(*discs)
        assert_equal(testmat2.todense(), desired)
        self.assert_(testmat2 is not testmat)
        self.assert_(isinstance(
            testmat, DiscretiserMatrices.DiscretiserMatrices.sparseType))        

    def test_mass(self):
        inst = DiscretiserMatrices.DiscretiserMatrices(self.disc1)
        desired_sparse = Matrices.sparseType(
            SystemMatrix.self_projection_matrix(self.disc1))
        desired = desired_sparse.todense()
        self._test_matrix(inst, desired, 'mass')
        from random import random
        # sparse lu objects can't be compared, so solve for a random vector
        vec = N.array([random() for x in xrange(desired.shape[0])])
        inst.mass_LU().solve(vec)
        assert_almost_equal(inst.mass_LU().solve(vec),
                            linsolve.factorized(desired_sparse)(vec),
                            decimal=14)
        
    def test_stiffness(self):
        desired = SystemMatrix.self_projection_matrix(self.disc1.D()).todense()
        inst = DiscretiserMatrices.DiscretiserMatrices(self.disc1)
        self._test_matrix(inst, desired, 'stiffness')

    def test_projectionOnto(self):
        desired = SystemMatrix.projection_matrix(
            self.disc2, self.disc1).todense()
        inst = DiscretiserMatrices.DiscretiserMatrices(self.disc1)
        self._test_matrix(inst, desired, 'projectionOnto', self.disc2)
        
    def test_exteriorDerivative(self):
        desired=N.zeros((7,9), N.int8)
        row = [1, 1, -1]
        desired[0,[0,3,1]] = row
        desired[1,[0,4,2]] = row
        desired[2,[1,5,2]] = row
        desired[3,[3,5,4]] = row
        desired[4,[1,7,6]] = row
        desired[5,[2,8,6]] = row
        desired[6,[5,8,7]] = row
        
        inst = DiscretiserMatrices.DiscretiserMatrices(self.disc1)
        self._test_matrix(inst, desired, 'exteriorDerivative', self.disc2)
        
