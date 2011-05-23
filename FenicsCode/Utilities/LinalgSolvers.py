from checkbox.job import PASS
__author__ = "Evan Lezar"
__date__ = "11 May 2011"


import numpy as np
import scipy.sparse
import scipy.sparse.linalg
from time import time

class SystemSolverBase ( object ):
    """
    A base class for the implementation of various solvers for sparse eigen systems.
    This base class provides logging functionality, but requires an extention to allow for actual solver implementation.
    """
    def _timestamp (self, id, res=np.nan ):
        """
        Add a timestamp/id pair to the logging data used to measure the progress of the solver
        """
        self._logging_data['time'].append ( time() )
        self._logging_data['id'].append ( id )
        self._logging_data['res'].append ( res )
    
    def set_preconditioner (self, M_type ):
        self._timestamp('preconditioner::start')
        if M_type is None:
            M = None
        elif M_type.lower() == 'diagonal':
            M = scipy.sparse.spdiags(1./self._A.diagonal(), 0, self._A.shape[0], self._A.shape[1])
        elif M_type.lower() == 'ilu':
            self._M_data = scipy.sparse.linalg.spilu ( self._A.tocsc() )
            M = scipy.sparse.linalg.LinearOperator ( self._A.shape, self._M_data.solve )
        else:
            print "Warning: Preconditioner type '%s' not recognised. Not using preconditioner." % M_type
            M = None 
        
        self._M = M
        self._timestamp('preconditioner::end')
        
    def __init__ ( self, A, preconditioner_type = None ):
        self._logging_data = { 'time': [], 'id': [], 'res': [] }
        self._timestamp( 'init' )
        self._A = A
        self.set_preconditioner ( preconditioner_type )
        self._callback_count = 0
        
    def _callback ( self, xk ):
        self._callback_count += 1;
        self._timestamp( self._callback_count, res=np.linalg.norm( self._A.matvec( xk ) - self._b ) )
    
    def set_b ( self, b ):
        print b.shape
        self._b = b
        if len(b.shape) > 1:
            if b.shape[0] > b.shape[1]:
                self._b = self._b[:,0]
            else:
                self._b = self._b[0,:]
        
        print self._b.shape
        
    def solve ( self, b ):
        self._timestamp( 'solve::begin' )
        self.set_b(b)
        
        x, info = self._call_solver ()

        if ( info > 0 ):
            "convergence to tolerance not achieved, in %d iterations" % info
        elif ( info < 0 ): 
            "illegal input or breakdown (%d)" % info
            
        self._timestamp( 'solve::end' )
        return x
    
    def get_logging_data (self):
        return self._logging_data
    
    def plot_convergence (self, x_is_time=False, show_plot=False, label=None, style='-'):
        import pylab as P
        
        y_data = np.log10( np.array(self._logging_data['res']) )
        if x_is_time:
            t0 = self._logging_data['time'][0]
            x_data = np.array(self._logging_data['time'], dtype=np.float64) - t0
        else:
            index = np.where(np.isfinite(y_data))[0]
            y_data = y_data[index]
            x_data = np.zeros ( index.shape )
            for i in range(len(index)):
                x_data[i] = self._logging_data['id'][index[i]]
        
        P.plot ( x_data, y_data, style, label=label )
        
        if show_plot:
            P.legend(loc='upper right')
            P.show()
    
    def _call_solver (self):
        raise Exception ( "solver driver not implemented")

class BiCGStabSolver ( SystemSolverBase ):
    """
    This solver implements an iterative Stabilised BICG solver using scipy
    """
    def _call_solver (self):
        """
        Solves the linear system (self._A)x = self._b
        
        @see: The SystemSolverBase class
        """
        return scipy.sparse.linalg.bicgstab(self._A, self._b, M=self._M, callback=self._callback )

class GMRESSolver ( SystemSolverBase ):
    """
    This solver implements an iterative GMRES solver using scipy
    """
    def _call_solver (self):
        """
        Solves the linear system (self._A)x = self._b
        
        @see: The SystemSolverBase class
        """
        return scipy.sparse.linalg.gmres(self._A, self._b, M=self._M, callback=self._callback )


def solve_sparse_system ( A, b ):
    """
    This function solves the sparse linear system Ax = b for A a scipy sparse matrix, and b a numpy or scipy array
    
    An incomplete LU preconditioner is used with the bicgstab iterative solver
    
    @param A: a square matrix 
    @param b: the RHS vector
    """
    
    solver = BiCGStabSolver ( A, 'ilu' )
    x = solver.solve(b)
    return x
    
