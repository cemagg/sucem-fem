from __future__ import division

from numpy.testing import NumpyTestCase, assert_array_equal, assert_array_almost_equal,\
     assert_almost_equal, assert_equal
import numpy as N
from NewCode import Utilities

class test_Struct(NumpyTestCase):
    def setUp(self):
        self.data=dict(attr1='lala', attr2='lala2')
        self.struct=Utilities.Struct(**self.data)

    def test_init(self):
        assert_equal((self.struct.attr1, self.struct.attr2), ('lala', 'lala2'))
        assert_equal(self.struct, self.data)
        struct = Utilities.Struct(self.data)
        assert_equal((struct.attr1, self.struct.attr2), ('lala', 'lala2'))
        assert_equal(struct, self.data)

    def test_setattr(self):
        self.struct.attr3 = 'abc'
        assert_equal(self.struct.attr3, 'abc')
        self.data['attr3'] = 'abc'
        assert_equal(self.struct, self.data)
        def tmp():
            self.struct.__getitem__ = 123
        self.assertRaises(AttributeError, tmp)
        def tmp():
            self.struct['__getitem__']=1
        self.assertRaises(KeyError, tmp)
        self.data['attr4'] = 'def'
        self.struct['attr4'] = 'def'
        assert_equal(self.struct, self.data)
        assert_equal(self.struct.attr4, 'def')
        
class test_CacheLast(NumpyTestCase):
    class TestClass(object):
        __metaclass__ = Utilities.CacheLast
        @Utilities.CacheLast.CachedMethod
        @Utilities.add_attr('newattr', 'I am an extra function attribute')
        def cached(self, x):
            return set([x])
        def nocache(self,x):
            return set([x])

    def test(self):
        tc1 = self.TestClass()
        tc2 = self.TestClass()
        assert_equal(tc1.cached.newattr, 'I am an extra function attribute')
        c1 = tc1.cached(1)
        assert_equal(c1, set([1]))
        self.assert_(c1 is tc1.cached(1))
        c1_2 = tc2.cached(1)
        self.assert_(c1 is not c1_2)    # They should be seperately cached
        assert_equal(tc2.cached(1), set([1]))
        tc1.cached.clearCache()
        self.assert_(c1 is not tc1.cached(1))
        self.assert_(c1_2 is tc2.cached(1)) # Clearing tc1's cache doesn't affect tc2
        c2 = tc1.cached(2)              # Test that new parms cause the cache to be flushed
        assert_equal(c2, set([2]))
        self.assert_(c2 is tc1.cached(2))
        assert_equal(tc2.cached(1), set([1]))

class test_rechain(NumpyTestCase):
    def test_rechain(self):
        l1 = [0,1,2]
        l2 = ['a', 'b', 'c']
        c = Utilities.rechain(l1, l2)
        # Same result as itertools.chain, except a fresh iterator is generated
        # whenever the last one has run out.
        assert_equal([x for x in c], l1 + l2)
        # Standard itertools.chain would generate an empty list here
        assert_equal([x for x in c], l1 + l2)

class test_gen_lagrange_polys(NumpyTestCase):
    def test_fns(self):
        interp_pts = [0, 1/3, 1/2, 1]
        fns = Utilities.gen_lagrange_polys(interp_pts)
        assert_array_almost_equal([fn(p) for fn, p in zip(fns, interp_pts)],
                                  [1. for p in interp_pts], decimal=15)
        assert_array_almost_equal([[fn(pz) for pz in interp_pts if pz != pi]
                                   for fn, pi in zip(fns, interp_pts)],
                                  N.zeros((len(interp_pts), len(interp_pts) - 1)),
                                  decimal=14)

class test_on_box_surf(NumpyTestCase):
    eps_g = 1e-9
    a,b,c = 1,2,3
    tnds_a = N.array([[0,0,.5], [0,.5,0], [0.5,0,0]]) # False
    tnds_b = N.array([[0,0,.5], [0,0,0], [0,.5,0]]) # True
    tnds_c = N.array([[0,0,0], [0,b,0], [0,b,c]]) # True
    tnds_d = N.array([[.5, .8, c], [a,b,c], [.7, 1.2, c]]) # True
    tnds_e = N.array([[.5, .8, 0], [a,b,0], [.7, 1.2, 0]]) # True
    tnds_f = N.array([[.5, .8, 0], [a,b,c], [.7, 1.2, 0]]) # False
    tnds_g = N.array([[.5, .8, 0], [a,b+.1,0], [.7, 1.2, 0]]) # False
    tnds_h = N.array([[.5, b, .8], [a,b,c], [.7, b, 1.2]]) # True
    tnds_i = N.array([[.5, 0, .8], [a,0,c], [.7, 0, 1.2]]) # True
    tnds_j = N.array([[a,0,.5], [a,0,0], [a,.5,0]]) # True
    tnds_k = N.array([[a,0,.5], [a,0,0], [a,.5,c+0.1]]) # False
    tnds_l = N.array([[a,0,.5], [a,-0.1,0], [a,.5,c]]) # False
    tnds_m = N.array([[a,0,.5], [a,b+0.1,0], [a,.5,c]]) # False
    tnds_n = N.array([[a,0,.5], [a,0,0], [a,.5,-0.1]]) # False
    def test_box(self):
        tp = Utilities.on_box_surf([0,0,0], [self.a,self.b,self.c], self.eps_g)
        self.assert_(tp(self.tnds_a) == False)
        self.assert_(tp(self.tnds_b) == True )
        self.assert_(tp(self.tnds_c) == True )
        self.assert_(tp(self.tnds_d) == True )
        self.assert_(tp(self.tnds_e) == True )
        self.assert_(tp(self.tnds_f) == False)
        self.assert_(tp(self.tnds_g) == False)
        self.assert_(tp(self.tnds_h) == True )
        self.assert_(tp(self.tnds_i) == True )
        self.assert_(tp(self.tnds_j) == True )
        self.assert_(tp(self.tnds_k) == False)
        self.assert_(tp(self.tnds_l) == False)
        self.assert_(tp(self.tnds_m) == False)
        self.assert_(tp(self.tnds_n) == False)
