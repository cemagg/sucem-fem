import numpy as N
from scipy import sparse
from scipy.sparse.linalg import dsolve as linsolve

from NewCode.Utilities import CacheLast, Struct
from NewCode.MatrixUtils import Matrices
from NewCode import SystemMatrix
from NewCode.DifferentialForm import Operations

def ZeroDValueError(fun):
    def newFun(*name, **kwargs):
        try: return fun(*name, **kwargs)
        except AttributeError, e:
            if e.args[0].find("'D'") != -1:
                raise ValueError, \
                      "The exterior derivative of an exterior derivative is per defn 0"
            else: raise
    return newFun

class DiscretiserMatrices(Matrices):
    
    def __init__(self, disc):
        self.disc = disc

    @CacheLast.CachedMethod
    def mass(self):
        print "Calculating mass matrix"
        return self._finalizeMatrix(
            SystemMatrix.projection_matrix(self.disc))
    
    @CacheLast.CachedMethod
    def mass_LU(self):
        return Struct(solve=linsolve.factorized(self.mass()))
    
    @CacheLast.CachedMethod
    @ZeroDValueError
    def stiffness(self):
        print "Calculating stiffness matrix"
        return self._finalizeMatrix(
            SystemMatrix.projection_matrix(self.disc.D()))

    @CacheLast.CachedMethod
    def boundarySurfaceMatrix(self):
        print "Calculating boundarySurfaceMatrix"
        return self._finalizeMatrix(
            SystemMatrix.boundary_matrix(self.disc))

    @CacheLast.CachedMethod
    def projectionOnto(self, target_disc):
        """
        Projection matrix of this discretiser onto target_disc

        The matrix P = DiscretiserMatricesInstance.projectionOnto(target_disc)
        represents the basis functions of self.disc tested by the basis
        functions of target_disc. This also means that P's row-index relates to
        target_disc and its column index relates to self.disc.

        """
        print "Calculating Projection Matrix"
        return self._finalizeMatrix(
            SystemMatrix.projection_matrix(target_disc, self.disc))

    @CacheLast.CachedMethod
    def lumpy_projectionOnto(self, target_disc):
        """
        Projection matrix of this discretiser onto target_disc

        The matrix P = DiscretiserMatricesInstance.projectionOnto(target_disc)
        represents the basis functions of self.disc tested by the basis
        functions of target_disc. This also means that P's row-index relates to
        target_disc and its column index relates to self.disc. This assumes
        that target_disc has an alt_rule integration rule set up that
        corresponds to self's

        """
        print "Calculating Projection Matrix"
        return self._finalizeMatrix(
            SystemMatrix.lumpy_projection_matrix(target_disc, self.disc))


    @CacheLast.CachedMethod
    @ZeroDValueError
    def exteriorDerivative(self, target_disc):
        """
        Calculate exterior div of self i.t.o. target n+1 form target_disc

        Note target_disc must be able to represent the ext diff exactly
        """
        return Operations.ext_diff_mat(
            self.disc, target_disc)

    @CacheLast.CachedMethod
    @ZeroDValueError
    def partialExteriorDerivative(self, target_disc):
        """
        Calculate exterior div of self i.t.o. target n+1 form target_disc

        target_disc only needs to be able to represent a subset of the ext diff
        exactly.
        """
        return Operations.ext_diff_mat(
            self.disc, target_disc, ignore_missing_2form=True)
