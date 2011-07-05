from __future__ import division

import numpy as np

def normalised_RMS(actual, desired, extra_error=None):
    """Calculate normalised RMS error %

    Parameters
    ----------

    actual -- actual values
    desired -- desired values, the series used for normalisation
    extra_error -- An optional extra orthogonal error component

    Return value

    The normalised RMS error of actual vs desired in percentage
    """
    err = np.abs(actual-desired)
    if extra_error is None:
        extra_error = np.zeros_like(err)
    else:
        extra_error = np.abs(extra_error)

    return 100*np.sum((err/np.abs(actual))**2 + (extra_error/np.abs(actual))**2)
     
        
