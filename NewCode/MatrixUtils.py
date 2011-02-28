from __future__ import division

import numpy as N
from scipy import sparse
from scipy.sparse.linalg import iterative
from scipy.sparse.linalg import dsolve as linsolve

from NewCode.Utilities import CacheLast

class IterativeSolverError(ArithmeticError):
    def __str__(self):
        return 'Illegal input or breakdown while solving.'
    
class ConvergenceError(IterativeSolverError):
    def __init__(self, no_iters):
        self.no_iters = no_iters

    def __str__(self):
        return "Solver failed to converge after %s iterations." % self.no_iters

class MatrixSolver(object):
    iterativeSolveTol=1e-9
    maxIterFac=1000
    def __init__(self, tol=None, *names, **kwargs):
        if tol: self.iterativeSolveTol = tol
        super(MatrixSolver, self).__init__(*names, **kwargs)

    def solve_mat_vec(self, A, b):
        """
        Solve matrix equation [A][x] = [b] for [x]

        This method wraps matrix solution calls so that different matrix
        solution strategies don't require wide-ranging code changes.
        """

        #print 'Solving Matrix'
        try:
            return A.solve(b)
        except AttributeError:
            pass
        (x, info) = iterative.cg(
            A,b,tol=self.iterativeSolveTol, maxiter=len(b)*self.maxIterFac)
        if info > 0:
            raise ConvergenceError, info
        if info < 0:
            raise IterativeSolverError
        return x

    def solve_mat_matvec(self, A, B, b):
        """
        Solve matrix equation [A][x] = [B][b] for [x]

        This method wraps matrix solution calls so that different matrix
        solution strategies don't require wide-ranging code changes.
        """

        return self.solve_mat_vec(A,B*(b))


def add_diagonal_preconditioner(A):
    M_inv = get_diag(A)
    def psolve(b):
        #print "psolve"
        return M_inv*b
    A.psolve = psolve

get_diag = lambda A: N.array([1/A[i,i][0] for i in xrange(A.shape[0])], A.dtype)


class Matrices(object):
    __metaclass__ = CacheLast
    sparseType = sparse.csc_matrix
    maintainSorted = True
    dtype = None

    def clearAll(self):
        for ca in self.cachedAttrs.values():
            ca.clearCache()

    def _finalizeMatrix(self, matr):
        if sparse.isspmatrix(matr):
            if self.dtype is None:
                matr = self.sparseType(matr)
            else:
                matr = self.sparseType(matr).astype(self.dtype)
            if self.maintainSorted:
                if not matr.has_sorted_indices: matr.sort_indices()
        return matr

def merge_sparse_blocks(block_mats, format='coo', dtype=N.float64):
    """
    Merge several sparse matrix blocks into a single sparse matrix

    Input Params
    ============

    block_mats -- sequence of block matrix offsets and block matrices such that
                  block_mats[i] == ((row_offset, col_offset), block_matrix)
    format -- Desired sparse format of output matrix
    dtype -- Desired dtype, defaults to N.float64

    Output
    ======

    Global matrix containing the input blocks at the desired block
    locations. If csr or csc matrices are requested it is ensured that their
    indices are sorted.

    Example
    =======

    The 5x5 matrix A containing a 3x3 upper diagonal block, A_aa and 2x2 lower
    diagonal block A_bb:

    A = [A_aa  0   ]
        [0     A_bb]

    A = merge_sparse_blocks( ( ((0,0), A_aa), ((3,3), A_bb)) )

    """
    
    nnz = sum(m.nnz for o,m in block_mats)
    data = N.empty(nnz, dtype=dtype)
    row = N.empty(nnz, dtype=N.intc)
    col = N.empty(nnz, dtype=N.intc)
    nnz_o = 0
    for (row_o, col_o), bm in block_mats:
        bm = bm.tocoo()
        data[nnz_o:nnz_o+bm.nnz] = bm.data
        row[nnz_o:nnz_o+bm.nnz] = bm.row+row_o
        col[nnz_o:nnz_o+bm.nnz] = bm.col+col_o
        nnz_o += bm.nnz

    merged_mat = sparse.coo_matrix((data, (row, col)), dtype=dtype)

    if format != 'coo':
        merged_mat = getattr(merged_mat, 'to'+format)()
        if format == 'csc' or format == 'csr':
            if not merged_mat.has_sorted_indices: merged_mat.sort_indices()

    return merged_mat

def extract_diag_mat(mat):
    return sparse.dia_matrix(([mat.diagonal()], [0]), shape=mat.shape, dtype=mat.dtype)

class EmptySparse(N.matrix):
    def __new__(subtype, shape=(0,0)):
        return N.asmatrix(N.zeros(shape))
    nnz = 0

