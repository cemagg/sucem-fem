__author__ = "Evan Lezar"
__date__ = "11 May 2011"


import numpy as np

def solve_sparse_system ( A, b ):
    """
    This function solves the sparse linear system Ax = b for A a scipy sparse matrix, and b a numpy or scipy array
    
    A inverse diagonal preconditioner is used with the bicgstab iterative solver
    
    @param A: a square matrix 
    @param b: the RHS vector
    """
    import scipy.sparse
    import scipy.sparse.linalg
    
    x, info = scipy.sparse.linalg.bicgstab( A, b, M=scipy.sparse.spdiags(1./A.diagonal(), 0, A.shape[0], A.shape[1]) )
    
    if ( info > 0 ):
        "convergence to tolerance not achieved, in %d iterations" % info
    elif ( info < 0 ): 
        "illegal input or breakdown (%d)" % info
    
    return x
    
