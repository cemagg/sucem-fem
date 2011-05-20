__author__ = "Evan Lezar"
__date__ = "11 May 2011"


import numpy as np





def solve_sparse_system ( A, b, calc_residuals=False, solver='bicgstab', preconditioner='ilu' ):
    """
    This function solves the sparse linear system Ax = b for A a scipy sparse matrix, and b a numpy or scipy array
    
    A inverse diagonal preconditioner is used with the bicgstab iterative solver
    
    @param A: a square matrix 
    @param b: the RHS vector
    """
    import scipy.sparse
    import scipy.sparse.linalg
    
    if calc_residuals:
        residual_list = []
        def myCallback ( xk ):
            residual_list.append( np.linalg.norm( A.matvec( xk ) - b[:,0] ) )
    
    if preconditioner == 'diagonal':
        M = scipy.sparse.spdiags(1./A.diagonal(), 0, A.shape[0], A.shape[1])
    elif preconditioner == 'ilu':
        M_tmp = scipy.sparse.linalg.spilu ( A.tocsc() )
        M = scipy.sparse.linalg.LinearOperator ( A.shape, M_tmp.solve )
    else:
        M = None;
    
    if solver == 'bicgstab':
        x, info = scipy.sparse.linalg.bicgstab(A, b, M=M, callback=myCallback )
    elif solver == 'gmres':
        x, info = scipy.sparse.linalg.gmres(A, b, M=M, callback=myCallback )
        
    if ( info > 0 ):
        "convergence to tolerance not achieved, in %d iterations" % info
    elif ( info < 0 ): 
        "illegal input or breakdown (%d)" % info
    
    if calc_residuals:
        return x, residual_list
    
    return x
    
