# Makes 1/2 (etc.) return floating point rather than 0.
from __future__ import division

from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal
from numpy import array, arange, float64, int32

from NewCode import Mesh
from NewCode.DifferentialForm.BasisFunction import Twoform
from NewCode.Utilities import Struct
from NewCode.tests.TestMeshes import FlatTet

class test_basisfunctions(TestCase):
    def setUp(self):
        self.local_facenodes = Mesh.Element.LOCAL_FACENODES
        
    def test_face0(self):
        #Yes, these don't add to 1, but they're fine for testing
        lam = array([1, 10, 100, 1000], float64)
        desired = array([[100,-10,0,1,0,0],
                         [1000,0,-10,0,1,0],
                         [0,1000,-100,0,0,1],
                         [0,0,0,1000,-100,10]], float64)*2
        assert_almost_equal([facefun(lam)
                             for facefun in
                             (facefun_g(self.local_facenodes)
                              for facefun_g in Twoform.facefuns0)],
                            desired, decimal=15)

    def test_face0_D(self):
        # lam is actually not used by the lowest order face_D funcs so just
        # pass 0 as a dummy value
        assert_almost_equal([facefun_D(0) for facefun_D in
                             (facefun_g(self.local_facenodes)
                              for facefun_g in Twoform.facefuns0_D)],
                            array([[6], [-6], [6], [-6]]), decimal=15)

class test_basis_set(TestCase):
    def test_1_mixed(self):
        bset = Twoform.basis_set(1, mixed=True)
        # Check default value of mixed
        assert_equal(Twoform.basis_set(1), bset)
        desired=Struct(fns=dict(face=Twoform.facefuns0),
                       fns_D=dict(face=Twoform.facefuns0_D),
                       info=Struct(form=2, type='mmbotha06_solh', order=0.5))
        assert_equal(bset, desired)
        self.assertRaises(AssertionError, Twoform.basis_set, 0)
        self.assertRaises(NotImplementedError, Twoform.basis_set, 1000000)


    def test_1(self):
        bset = Twoform.basis_set(1, mixed=False)
        desired=Struct(fns=dict(face=Twoform.facefuns0+Twoform.facefunsR1),
                       fns_D=dict(face=Twoform.facefuns0_D+Twoform.facefunsR1_D),
                       info=Struct(form=2, type='mmbotha06_solh', order=1))
        assert_equal(bset, desired)

    def test_2(self):
        bset = Twoform.basis_set(2,mixed=False)
        desired=Struct(fns=dict(face=Twoform.facefuns0+Twoform.facefunsR1
                                +Twoform.facefunsR2,
                                vol=Twoform.volfunsR2),
                       fns_D=dict(face=Twoform.facefuns0_D+Twoform.facefunsR1_D
                                  +Twoform.facefunsR2_D,
                                vol=Twoform.volfunsR2_D),
                       info=Struct(form=2, type='mmbotha06_solh', order=2))
        assert_equal(bset, desired)
