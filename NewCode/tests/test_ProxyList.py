from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal
from itertools import izip
import sys
#
# Local Imports
#
from NewCode import ProxyList

class test_ProxyList(NumpyTestCase):
    def setUp(self):
        self.test_attrs = {'colour' : ['blue', 'black', 'white', 'mauve'],
                      'id' : [1,2,3,4]}
        self.len = 4
        self.TestItemClass = ProxyList.ItemClassFactory(self.test_attrs,
                                                        'TestItemClass')
        self.testItemInstance = self.TestItemClass(self.test_attrs)
        self.test_list = ProxyList.ProxyList(self.testItemInstance)
        
    def check_counts(self, level=1):
        cnt = len([x for x in self.test_list])
        # Is the correct number if items iterated over?
        assert_equal(cnt, self.len)
        # Does the __len__ method return the correct length?
        assert_equal(len(self.test_list), self.len)
        # Are the correct bounds stored by the list class?
        assert_equal(self.test_list.bounds, self.len-1)

    def check_values(self, level=1):
        # Using iterators
        assert_equal([x.colour for x in self.test_list],
                     self.test_attrs['colour'])
        assert_equal([x.id for x in self.test_list],
                     self.test_attrs['id'])
        # Using subscripts
        for i in range(self.len):
            assert_equal(self.test_list[i].colour,
                         self.test_attrs['colour'][i])
            assert_equal(self.test_list[i].id,
                         self.test_attrs['id'][i])

    def check_slices(self, level=1):
        assert_equal([x.colour for x in self.test_list[1:3]],
                      self.test_attrs['colour'][1:3])

    def test_slice_len(self):
        assert_equal(len(self.test_list[1:100]), 3)
        
    def check_attr(self):
        proxy_attrs = tuple(self.test_attrs.keys())
        assert_equal(self.test_list.entity.proxy_attrs,
                     proxy_attrs)

    def check_list_repr(self):
        assert_equal(self.test_list.list_repr(),
                     self.test_attrs)


class test_Memoize_index(NumpyTestCase):
    def setUp(self):
        self.test_attrs = {'colour' : ['blue', 'black', 'white', 'mauve'],
                           'id' : [1,2,3,4]}
        self.len = 4
        self.TestItemBaseClass = ProxyList.ItemClassFactory(self.test_attrs,
                                                            'TestItemBaseClass')
        TestItemBaseClass = self.TestItemBaseClass
        class TestItemClass(TestItemBaseClass):
            def __init__(self, *names, **kwargs):
                self._lastindex = None
                super(TestItemClass,self).__init__(*names, **kwargs)

            @ProxyList.memoize_last(args=None)
            def id_squared(self):
                assert(self.index != self._lastindex)
                self._lastindex = self.index
                return self.id**2

        self.TestItemClass = TestItemClass
        self.testItemInstance = self.TestItemClass(self.test_attrs)
        self.test_list = ProxyList.ProxyList(self.testItemInstance)


    def test_Memoize_index(self):
        desired = [1,4,9,16]
        assert_equal([(item.id_squared(), item.id_squared())
                      for item in self.test_list],
                     zip(desired, desired))

    def test_Memoize_multiple_instances(self):
        desired1 = [1,4,9,16]
        desired2 = [16,4,1,9]
        test_list1 = self.test_list
        testItemInstance2 = self.TestItemClass({'colour': self.test_attrs['colour'],
                                                'id': [4,2,1,3]})
        test_list2 = ProxyList.ProxyList(testItemInstance2)
        assert_equal([(id1.id_squared(), id2.id_squared())
                      for id1, id2 in izip(test_list1, test_list2)],
                     zip(desired1, desired2))

