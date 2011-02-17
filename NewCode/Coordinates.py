from __future__ import division
"""
Mix-in Classes for dealing with coordinate systems. 

Module Classes
==============

SimplexCoord
  Deals with tetrahedral simplex coordinates

"""

import numpy as N
from numpy import array, float64, ones, zeros, linalg, arange, compress, cross, \
     mod, newaxis, dot
from NewCode import ProxyList


class SimplexCoord(object):
    noLocalCoords = 4
    """
    Mix-in class for simplex tetrahedral coordinates.

    Mix-in Interface
    ================

    Attributes used
    ---------------

    self.nodeCoords
      4-row x 3-col array, each row containing the x,y,z coords of a node. If
      the mix-in class does not provide this, it can also be passed to the
      __init__ method

    self.LOCAL_FACENODES
      4-row x 3-col array, each row containing the three local node indices
      that define the face.

    Methods
    =======

    __init__
      The element nodeCoords attribute can be specified if the mix-in class lacks
      it.

    GradLambda
      Calculate the gradient of the four simplex coordinates.
      
    covBaseVecs
      Calculate the Covariant unitary vectors. Aliased to GradLambda, since
      they are the same for simplex tets
    
    FaceBasisVecs
      Calculate the 2-form basis vectors for a given face.
    """

    def __init__(self, nodeCoords=None, *names, **kwargs):
        super(SimplexCoord, self).__init__(*names, **kwargs)
        if nodeCoords is not None: self.nodeCoords = nodeCoords
        try:                            # for @ProxyList.method compatibility
            getattr(self, 'index')
        except AttributeError:
            self.index = 0

    def simplexMat(self):
        simplexMat = ones((4,4), float64)
        simplexMat[0:4, 1:4] = self.nodeCoords
        return simplexMat

    @ProxyList.memoize_last(args=None)
    def GradLambda(self):
        """
        Calculate the gradient of the four simplex coordinates

        Arguments:

        elno -- The element number
        geomwrap -- f2py wrapped FEMFEKO geometry data module

        Returns:

        4x3 matrix, with each row containing the x, y and z components of one
        simplex coordiate. I.e. grad_lamda_tet[i, :] are the components for
        coordinate i.

        """
        #
        # The gradients are calculated as follows:
        #
        #                             ^  ^  ^ 
        #   __               +1   | 0 x  y  z  | 
        #   \/ lambda_1  =   --   | 1 x2 y2 z2 |
        #                    6V   | 1 x3 y3 z3 | 
        #                         | 1 x4 y4 z4 |
        #
        # and similar expressions for the others (the signs alternate):
        #                             ^  ^  ^ 
        #   __               -1   | 0 x  y  z  | 
        #   \/ lambda_2  =   --   | 1 x1 y1 z1 |
        #                    6V   | 1 x3 y3 z3 | 
        #                         | 1 x4 y4 z4 |
        # 
        # Note there that the volume(V) here may be SIGNED!
        #
        # We do this i.t.o the simplex coordinate matrix:
        #
        #         [ 1 x1 y1 z1 ] 
        #         [ 1 x2 y2 z2 ] 
        #         [ 1 x3 y3 z3 ] 
        #         [ 1 x4 y4 z4 ]
        #
        #         where x1 refers to the x-coord of node 1 of the tet, etc.

        simplexMat = self.simplexMat()

        #
        # volFactor = 1/6V, with 6V = det(simplexMat)
        volFactor = 1/linalg.det(simplexMat)
        #
        # Initialise the gradients to zero
        grads = zeros((4,3), float64)

        # Now we will use co-factor expansion accross the top row of the matrices
        # described above to calculate the determinants. This let's us calculate
        # each Cartesian component individually.
        #
        # The co-factor matrices can be constructed by removing the row of
        # simplexMat corresponding to simplex co-ord being calculated, and the
        # column corresponding to the cartesian co-ordinate being calculated
        #
        # E.g. when calculating the x-compoment of del(lambda_1), remove the row
        # containing [1 x1 y1 z1], and the column containing [ x2 x3 x4 ]

        ind = arange(4)
        sign = -1
        for var in range(4):                # Loop over the 4 simplex variables
            for component in range(3):      # Loop over x, y, z
                row = var                   
                col = component+1
                detmat = compress(ind != row,
                                  compress(ind != col, simplexMat, axis=1),
                                  axis=0)
                grads[var, component] = linalg.det(detmat)*sign
                sign *= -1                  # Alternate the sign

        grads *= volFactor                  # Multiply by 1/6V
        return grads

    def covBaseVecs(self):
        """
        Basis vectors for Cartesian representation of covariant components.

        These are the 'Reciprocal Unitary Vectors' accoring to Stratton41's notation
        """
        return self.GradLambda().transpose()

    @ProxyList.memoize_last(args='one')
    def FaceBasisVecs(self,local_faceno):
        """
        Calculate the 2-form basis vectors for a given face.

        Arguments
        =========

        local_faceno
          The element local face number

        Output
        ======

        3 x 3 array, with each row containing the x, y, z component of a face
        basis vector.

        For 2-form face functions as defined in Lee97, the three basis vectors
        depend on the nodes that define the face. For a face defined by nodes
        A, B and C, the vectors are

        grad(lambda_B) x grad(lambda_C), grad(lambda_C) x grad(lambda_A),
        grad(lambda_A) x grad(lambda_B)

        Nodes A, B and C are determined according to the local face-node
        numbering as stored in self.LOCAL_FACENODES

        """
        cov = self.covBaseVecs().transpose()
        conv = zeros((3,3), float64)
        facenodes = self.LOCAL_FACENODES[local_faceno]
        for i in range(3):
            conv[i,:] = cross(cov[facenodes[mod(i+1, 3)], :],
                              cov[facenodes[mod(i+2, 3)], :])

        return conv

    def conBaseVecs(self):
        """
        Calculate matrix of all possible contravariant component basis vectors

        Not scaled by the magnitude of the transform's Jacobian determinant.

        The matrix is 3x6, where ConUnitVecs[:,k] = grad(lambda_i) x grad(lamda_j)
        for (i,j) in the sequence ((0,1), (0,2), (0,3), (1,2), (1,3), (2,3))

        These are 'Unitary Vectors' accoring to Stratton41's notation.
        """
        cov = self.covBaseVecs()
        return array([cross(cov[:, i], cov[:, j])
                       for i,j in ((0,1), (0,2), (0,3), (1,2), (1,3), (2,3))],
                     float64).transpose()
    
    def local2global(self,lam):
        """
        Return the global coordinates of a given local coordinate point
        """
        return sum(self.nodeCoords*lam[:, newaxis]) 
            
    def global2local(self, xyz):
        """
        Return the local coordinates of a given global coordinate point
        """
        simplexMat = self.simplexMat()
        # Use the transpose-inverse of simplex matrix to calculate the simplex
        # coordinates from xyz. By Cramer's rule this is equivalent to
        # calculating the volume-ratio determinant of each simplex coord.
        xyz1 = [1]
        xyz1.extend(xyz)
        return dot(linalg.inv(simplexMat.T), xyz1)

    def face_coords2vol_coords(self, faceno, faceCoords):
        o = self.LOCAL_FACE_OPPOSING_NODE[faceno]
        lfc = list(faceCoords)
        lfc.insert(o, 0)
        return array(lfc, faceCoords.dtype)

    def vol_coords2face_coords(self, faceno, volCoords):
        o = self.LOCAL_FACE_OPPOSING_NODE[faceno]
        lvc = list(volCoords)
        del(lvc[o])
        return array(lvc, volCoords.dtype)
        
    def J(self):
        """
        Return the elemental Jacobian matrix

        For simplex coordinates defined as:

        [(r1-r4) (r2-r4) (r3-r4)]

        where r1 thru r4 are the position column vectors of the four tet-vertices
        """
        # TESTME!
        tmp = self.nodeCoords.T
        return tmp[:,0:3] - tmp[:,3,newaxis]

    def J_face(self, face_no):
        """
        Return the face Jacobian matrix

        For face with nodes a,b,c on tet with fourth node o,defined as:

        [(ra-ro) (rb-ro) (rc-ro)]

        where ra,rb,...  are position column vectors.
        """
        # TESTME!
        tmp = self.nodeCoords.T
        return tmp[:, self.LOCAL_FACENODES[face_no]] \
               - tmp[:,self.LOCAL_FACE_OPPOSING_NODE[face_no],newaxis]
        
    def InElement(self, xyz):
        """
        Returns True if the given point is in the element, False
        otherwise. Note that due to numerical tollerace issues it will return
        true also for points on the element boundary, and perhaps even a little
        outside of the element. This means that the same point can be 'found'
        in multiple elements. It is up to the user of this routine to decide
        how to deal with this situation.
        """

        # this: return abs(abs(self.global2local(xyz)).sum() - 1.0) < 4e-16 is
        # the eMAGUS style of determining if the point is in the element. It
        # requires a somewhat larger tollerance factor, though wether it is
        # truely more or less susceptible to numerical tollerances is hard to
        # say...

        # If all the local coordinates are positive the point is inside the
        # element. The ugly number relates to machine precision. Note, this is
        # _nasty_.

        return (self.global2local(xyz) >= -5e-16 ).all()

    
        
class BrickCoord(object):
    noLocalCoords = 6
    def __init__(self, nodeCoords=None, *names, **kwargs):
        super(BrickCoord, self).__init__(*names, **kwargs)
        if nodeCoords is not None: self.nodeCoords = nodeCoords
        try:                            # for @ProxyList.method compatibility
            getattr(self, 'index')
        except AttributeError:
            self.index = 0
    def covBaseVecs(self):
        return N.eye(3, dtype=N.float64) / self.gridStepSize
    def conBaseVecs(self):
        return N.eye(3, dtype=N.float64) / (
            self.gridStepSize[[1,2,0]]*self.gridStepSize[[2,0,1]])
    def J(self):
        return N.eye(3, dtype=N.float64) * self.gridStepSize

    def local2global(self, local_coords):
        # Good to test, but slow
        # assert(N.allclose(local_coords[0:3], 1-local_coords[3:6]))
        # Ignore the dependent coordinates by doing local_coords[0:3]
        nc = self.nodeCoords
        # Single or multiple elements
        n = nc[0] if len (nc.shape) == 2 else nc[:,0] 
        return n + local_coords[0:3]*self.gridStepSize
    
    def global2local(self, global_coord):
        half_coords = (global_coord - self.nodeCoords[0])/self.gridStepSize
        return N.array([half_coords, 1-half_coords],
                       half_coords.dtype).reshape((6,))

    def face_coords2vol_coords(self, faceno, faceCoords):
        insert = (0.,1.)
        lfc = faceCoords.tolist()
        lfc.insert(faceno % 3, insert[faceno//3])
        half_coords = N.array(lfc, faceCoords.dtype)
        return N.array([half_coords, 1-half_coords],
                       half_coords.dtype).reshape((6,))

    def vol_coords2face_coords(self, faceno, volCoords):
        lvc = volCoords[0:3].tolist()
        del(lvc[faceno % 3])
        return N.array(lvc, volCoords.dtype)

class PyramCoord(object):
    noLocalCoords = 5
    def covBaseVecs(self):
        nc = self.nodeCoords
        l1, l2, l3 = l123 = nc[[2,1,4]] - nc[0]
        J = N.dot(l1, N.cross(l2,l3))
        return N.array([N.cross(l123[(i+1)%3], l123[(i+2)%3])/J
                        for i in range(3)],
                       dtype=nc.dtype).T
         

    def conBaseVecs(self):
        nc = self.nodeCoords
        l1, l2, l3 = l123 = nc[[2,1,4]] - nc[0]
        J = N.dot(l1, N.cross(l2,l3))
        return l123.T/J

    def local2global(self, local_coords):
        ### TESTME!!!
        n0, n1, n2, n3, n4 = self.nodeCoords
        lam1,lam2,lam3,lam4,lam5 = local_coords
        return (lam1*lam2/(1-lam5)*n3 + lam1*lam4/(1-lam5)*n2
                + lam3*lam4/(1-lam5)*n0 + lam2*lam3/(1-lam5)*n1
                + lam5*n4)

        
