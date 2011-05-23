__author__ = "Evan Lezar"
__date__ = "11 May 2011"

"""this is a set of test cases for the testing of the linear algebra solvers in FenicsCode/Utilities"""

import sys
import unittest
import numpy as np


sys.path.insert(0, '../')
from FenicsCode.Utilities.LinalgSolvers import solve_sparse_system
del sys.path[0]


class TestSparseSolver ( unittest.TestCase ):
    def test_sparse_identity ( self ):
        import scipy.sparse
        N = 1000;
        A = scipy.sparse.eye ( N, N )
        b = np.random.rand ( N )
        
        x = solve_sparse_system ( A, b )
        
        np.testing.assert_array_equal( b, x )
        
        
if __name__ == "__main__":
    unittest.main()