__author__ = "Evan Lezar"
__date__ = "20 May 2011"

"""
This is a collection of routines to save and load finite element matrices to disk
"""
import os
import numpy as np

def check_path ( path, create=True ):
    """
    Check to see if a path exists, and create it if desired
    
    @param path: the path whose existence must be checked
    @param create: True by default. Create the path if it doesn't exist
    
    @return: True or False depending on whether the path exists after the function has executed
    """
    if not os.path.exists( path ):
        if not create:
            return True
        os.makedirs( path )
    return True

def save_scipy_matrix_as_mat ( path, name, matrix ):
    """
    Save a scipy sparse matrix as a .mat file.
    
    @param path: the folder in which the matrix is to be saved
    @param name: the filename of the matrix
    @param matrix: the scipy matrix to save
    """
    import scipy.io
    if not check_path ( path ):
        return False
    
    scipy.io.savemat ( os.path.join(path, name), {name: matrix }, oned_as='column' )
    
    return True


def load_scipy_matrix_from_mat ( path, name ):
    """
    Load a scipy sparse matrix from a .mat file.
    
    @param path: the folder in which the matrix is saved
    @param name: the filename of the matrix
    """
    import scipy.io
    import scipy.sparse
    
    filename =  os.path.join ( path, name)
    
    if not os.path.exists( filename ):
        return None
    
    data = scipy.io.loadmat ( filename )[name]
    
    if type(data) is np.ndarray:
        matrix = data
    else:
        matrix = data.tocsr ()
         
    return matrix
    
    
    