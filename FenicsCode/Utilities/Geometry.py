from __future__ import division
import numpy as np

def unit_vector(vec):
    """Calculate a unit vector of vec"""
    return vec/vector_length(vec)

def vector_length(vec):
    """Calculate the length of a vector"""
    return np.sqrt(np.sum(vec**2))
