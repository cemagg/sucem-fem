from __future__ import division
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal
from numpy import array
import numpy as N
from scipy import sparse
from scipy.sparse.linalg import iterative
from NewCode import MatrixUtils

class test_MatrixSolver(TestCase):

    def setUp(self):
        self.tol = 1e-10
        self.inst = MatrixUtils.MatrixSolver(self.tol)
        self.A = sparse.dok_matrix(array([[1.,2],[2.,1.]]))
        self.B = sparse.dok_matrix(array([[1.,2,3],[5,7,9]]))
        self.b = array([1.,2,3])

    def test_init(self):
        assert_equal(self.inst.iterativeSolveTol, self.tol)

    def test_solve_mat_vec(self):
        assert_array_equal(self.inst.solve_mat_vec(self.A, self.b[0:2]),
                           iterative.cg(self.A, self.b[0:2], tol=self.tol)[0])
        self.A[0,[0,1]]=0
        self.assertRaises(MatrixUtils.ConvergenceError,
                          self.inst.solve_mat_vec, self.A, self.b[0:2]
                          )

    def test_solve_mat_matvec(self):
        assert_equal(self.inst.solve_mat_matvec(self.A, self.B, self.b),
                     iterative.cg(self.A,
                                  self.B*(self.b),
                                  tol=self.tol)[0]
                     )
        
class test_merge_sparse_blocks(TestCase):
    A = N.array([[1  ,2,  0,2.5,1.1],
                 [0  ,3,1.5,  0,3.5],
                 [2.1,0,  4,  1,1.3],
                 [0,  0,  0,2.4,2.5],
                 [0,  0,  0,1.2,3.7]], N.float64)
    aa = sparse.lil_matrix(A[0:3, 0:3])
    ab = sparse.coo_matrix(A[0:3, 3:5])
    bb = sparse.csr_matrix(A[3:5, 3:5])

    def test_merge(self):
        merged_mat = MatrixUtils.merge_sparse_blocks((
            ((0,0), self.aa), ((0,3), self.ab), ((3,3), self.bb)))
        self.assert_(merged_mat.format == 'coo')
        assert_equal(merged_mat.todense(), self.A)

    def test_merge_to_csc(self):
        merged_mat = MatrixUtils.merge_sparse_blocks((
             ((3,3), self.bb), ((0,3), self.ab), ((0,0), self.aa)),
                                                     format='csc')
        self.assert_(merged_mat.format == 'csc')
        self.assert_(merged_mat.has_sorted_indices)
        assert_equal(merged_mat.todense(), self.A)

    def test_merge_to_csr(self):
        merged_mat = MatrixUtils.merge_sparse_blocks((
             ((3,3), self.bb), ((0,3), self.ab), ((0,0), self.aa)),
                                                     format='csr')
        self.assert_(merged_mat.format == 'csr')
        self.assert_(merged_mat.has_sorted_indices)
        assert_equal(merged_mat.todense(), self.A)
