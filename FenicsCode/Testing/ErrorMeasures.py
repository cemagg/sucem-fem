from __future__ import division

import numpy as np

def normalised_RMS(actual, desired, extra_error=None):
    """Calculate normalised RMS error %

    @param actual: actual values
    @param desired: desired values, the series used for normalisation
    @param extra_error: An optional extra orthogonal error component

    @return: The normalised RMS error of actual vs desired in percentage
    """
    err = np.abs(actual-desired)
    if extra_error is None:
        extra_error = np.zeros_like(err)
    else:
        extra_error = np.abs(extra_error)

    return 100*np.sum((err/np.abs(actual))**2 + (extra_error/np.abs(actual))**2)
     
        
def max_normalised_RMS(actual, desired, extra_error=None):
    """Calculate RMS error % normalised to the maximum value in actual

    This is to prevent values close to zero resulting in a large error
    norm while the absolute error remains small

    @param actual: actual values
    @param desired: desired values, the series used for normalisation
    @param extra_error: An optional extra orthogonal error component
    @return: The normalised (wrt. the maximum) RMS error of actual vs desired in percentage
    """
    err = np.abs(actual-desired)
    if extra_error is None:
        extra_error = np.zeros_like(err)
    else:
        extra_error = np.abs(extra_error)

    return 100*np.sum((err/np.abs(actual).max())**2 + 
                      (extra_error/np.abs(actual).max())**2)

