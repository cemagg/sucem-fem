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
    def __init__ ( self, A, preconditioner_type=None ):
        """
        The constructor for a System Solver
        
        @param A: The matrix for the system that must be solved
        @param preconditioner_type: A string indicating the type of preconditioner to be used
        """
        self._logging_data = { 'time': [], 'id': [], 'res': [] }
        self._timestamp( 'init' )
        self._A = A
        self.set_preconditioner ( preconditioner_type )
        self._callback_count = 0
    
    def _call_solver (self):
        """
        A stub routine that is to be implemented by an inherited class depending on the solver being used.
        """
        raise Exception ( "Solver driver not implemented.")
    
    def _callback ( self, xk ):
        """
        A simple callback routine used to track the progress of an iterative solver for Ax = b
        
        @param xk: the solution vector x at step k
        """
        self._callback_count += 1;
        if type(xk) is float:
            self._timestamp( self._callback_count, res=xk )
        else:
            self._timestamp( self._callback_count, res=calculate_residual ( self._A, xk, self._b ) )
        
    def _timestamp (self, id, res=np.nan ):
        """
        Add a timestamp/id pair to the logging data used to measure the progress of the solver
        
        @param id: an identifier for the timestamp
        @param res: NaN by Default. An optional residual parameter to store along with the timestamp info
        """
        self._logging_data['time'].append ( time() )
        self._logging_data['id'].append ( id )
        self._logging_data['res'].append ( res )
    
    def set_b ( self, b ):
        """
        Set the right-hand side vector b in the linear system Ax = b
        
        @param b: the right-hand side vector
        """
        self._b = b

    def set_preconditioner (self, M_type ):
        """
        Set the preconditioner (self._M) used in the solver
        
        @param M_type: A string to specify the preconditioner used
        """
        self._timestamp('preconditioner::start')
        if M_type is None:
            M = None
        elif M_type.lower() == 'diagonal':
            M = scipy.sparse.spdiags(1./self._A.diagonal(), 0, self._A.shape[0], self._A.shape[1])
        elif M_type.lower() == 'ilu':
            self._M_data = scipy.sparse.linalg.spilu ( self._A.tocsc(), drop_tol=1e-8, fill_factor=1  )
            M = scipy.sparse.linalg.LinearOperator ( self._A.shape, self._M_data.solve )
        else:
            print "Warning: Preconditioner type '%s' not recognised. Not using preconditioner." % M_type
            M = None 
        
        self._M = M
        self._timestamp('preconditioner::end')
        
        
    def solve ( self, b ):
        """
        Solve the linear system Ax = b for x
        
        @param b: the right-hand side vector
        """
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
        y_label = 'Residual [log10]'
        if x_is_time:
            t0 = self._logging_data['time'][0]
            x_data = np.array(self._logging_data['time'], dtype=np.float64) - t0
            x_label = 'Time [s]'
        else:
            index = np.where(np.isfinite(y_data))[0]
            y_data = y_data[index]
            x_data = np.zeros ( index.shape )
            for i in range(len(index)):
                x_data[i] = self._logging_data['id'][index[i]]
            x_label = 'Iterations'
            
        P.plot ( x_data, y_data, style, label=label )
        
        P.xlabel ( x_label )
        P.ylabel ( y_label )
        
        P.grid ( True )
        if show_plot:
            P.legend(loc='upper right')
            P.show()
    
    
    def _get_time_for_id (self, id ):
        return self._logging_data['time'][self._logging_data['id'].index(id)]
    
    def _get_elapsed_time ( self, id0, id1, force_total=False ):
        if force_total:
            return self._logging_data['time'][-1] - self._logging_data['time'][0]
         
        return self._get_time_for_id(id1) - self._get_time_for_id(id0)
        
    def print_logging_data ( self ):
        for k in self._logging_data:
            print k
            print self._logging_data[k]
    
    def print_timing_info ( self ):
        """
        Output the information contained in the logging data
        """
        
            
        solve_time = self._get_elapsed_time('solve::begin', 'solve::end')
        print 'solve time:', solve_time
        preconditioner_time = self._get_elapsed_time('preconditioner::start', 'preconditioner::end')
        print 'preconditioner time:', preconditioner_time
        
        print 'total time:', self._get_elapsed_time(None, None, True)
        
        
        

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

class UMFPACKSolver ( SystemSolverBase ):
    """
    This solver uses the UMFPACK fuctionality provided by scipy
    """
    def _call_solver (self):
        """
        Solves the linear system (self._A)x = self._b
        
        @see: The SystemSolverBase class
        """
        scipy.sparse.linalg.use_solver(useUmfpack=True)
        return scipy.sparse.linalg.spsolve ( self._A, self._b ), 0
    
    def plot_convergence (self, x_is_time=False, show_plot=False, label=None, style='-'):
        print "Direct solver has no convergence history"

class PyAMGSolver ( SystemSolverBase ):
    """
    This solver uses PyAMG to solve the linear system iteratively
    """
    def _call_solver (self):
        """
        Solves the linear system (self._A)x = self._b
        
        @see: The SystemSolverBase class
        """
        import pyamg
        x = pyamg.solve ( self._A, self._b, verb=True, tol=1e-8, maxiter=800 )
        return x, 0
    
    def plot_convergence (self, x_is_time=False, show_plot=False, label=None, style='-'):
        print "PAMG solver convergence display is not yet implemented"




def calculate_residual ( A, x, b ):
    """
    Calculate the residual of the system Ax = b
    
    @param A: a matrix
    @param x: a vector
    @param b: a vector 
    """
    return np.linalg.norm( A*x - b.reshape(x.shape) )


def solve_sparse_system ( A, b, preconditioner_type='ilu' ):
    """
    This function solves the sparse linear system Ax = b for A a scipy sparse matrix, and b a numpy or scipy array
    
    By default an incomplete LU preconditioner is used with the bicgstab iterative solver
    
    @param A: a square matrix 
    @param b: the RHS vector
    @param preconditioner_type: Preconditioner type string
    """
    
    solver = BiCGStabSolver ( A, preconditioner_type )
    x = solver.solve(b)
    return x
    
