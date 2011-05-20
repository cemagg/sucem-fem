__author__ = "Evan Lezar"
__date__ = "20 May 2011"

"""this is a set of test cases for the testing of the matrix IO routines in FenicsCode/Utilities"""

import numpy as np
import sys
import unittest

sys.path.insert(0, '../')
from FenicsCode.Utilities.MatrixIO import (
                                           load_scipy_matrix_from_mat,
                                           save_scipy_matrix_as_mat,
                                           )
del sys.path[0]

data_path = 'data/MatrixIO'


class TestMatrixIO ( unittest.TestCase ):
    def test_save_scipy_matrix ( self ):
        import scipy.sparse
        N = 1000;
        A = scipy.sparse.rand ( N, N, format='csr' )
        
        save_scipy_matrix_as_mat ( data_path, 'A', A )
    
    def __save_and_load_test (self, name, A ):
        save_scipy_matrix_as_mat ( data_path, name, A )
        A_load = load_scipy_matrix_from_mat ( data_path, name )
        
        np.testing.assert_equal ( A_load.shape, A.shape )
        
        if type(A) != np.ndarray:
            np.testing.assert_equal( A_load.data, A.data ) 
            np.testing.assert_equal( A_load.indices, A.indices )
            np.testing.assert_equal( A_load.indptr, A.indptr )
        else:
            np.testing.assert_equal( A_load, A )
        
            
        
    def test_save_and_load_real_matrix ( self ):
        import scipy.sparse
        N = 1000;
        A = scipy.sparse.rand ( N, N, format='csr' )
        self.__save_and_load_test( 'A_save_and_load_real', A)
        
    
    def test_save_and_load_complex_matrix ( self ):
        import scipy.sparse
        N = 1000;
        A = scipy.sparse.rand ( N, N, format='csr' ) + 1j*scipy.sparse.rand ( N, N, format='csr' )
        self.__save_and_load_test( 'A_save_and_load_complex', A)
    
    def test_save_and_load_1D_column ( self ):
        N = 1000;
        A = np.random.rand ( N,1 )
        self.__save_and_load_test( 'A_save_and_load_1D', A)
    
    def test_save_and_load_1D_row ( self ):
        N = 1000;
        A = np.random.rand ( 1,N )
        self.__save_and_load_test( 'A_save_and_load_1D', A )
                                   
        
if __name__ == "__main__":
    unittest.main()