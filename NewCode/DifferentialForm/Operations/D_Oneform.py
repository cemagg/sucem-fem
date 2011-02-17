import gc
from itertools import izip
import numpy as N
from scipy import sparse, linalg

from NewCode.Utilities import rechain
from NewCode.Exceptions import AllZero, DeadElement
from NewCode.MatrixUtils import Matrices                            

class D_oneformBase(object):
    def _test_basis_types(self, disc1, disc2):
        if disc1.basisSet.info.type != self.oneformType:
            raise NotImplementedError(
                'Wrong 1-form basis function type, '+self.oneformType+' expected')
        if disc2.basisSet.info.type != self.twoformType:
            raise NotImplementedError(
                'Wrong 2-form basis function type, '+self.twoformType+' expected')


    def _get_localMap(self, disc1, disc2):
        try: 
            return N.ascontiguousarray(self.curlLocalMaps[
                disc1.basisSet.info.order, disc2.basisSet.info.order].T)
        except KeyError:raise NotImplementedError(
            'Curl of %f order 1-form onto %f order 2-form not currently supported'
            % (disc1.basisSet.info.order, disc2.basisSet.info.order))

    def __call__(self, disc1, disc2, ignore_missing_2form=False):
        assert(disc1.mesh is disc2.mesh)
        assert(disc1.basisSet.info.form == 1)
        assert(disc2.basisSet.info.form == 2)
        self._test_basis_types(disc1, disc2)
        localMap = self._get_localMap(disc1, disc2)
        no_2fdofs, no_1fdofs = disc2.totalDOFs, disc1.totalDOFs
        mat = [dict() for i in range(no_2fdofs)]
        dofnos_1, dofnos_2, weights = localMap
        p1, p2 = disc1.permuter, disc2.permuter
        el_in_disc_1 = disc1.geomEntities.el_in_disc
        el_in_disc_2 = disc2.geomEntities.el_in_disc
        print 'Constructing Curl dict representation'
        for el1f, el2f in izip(disc1.elements, disc2.elements):
            if not el_in_disc_1(el1f) and not el_in_disc_2(el2f) : continue
            try:
                (l1,g1) = p1.permuteElement(el1f, remove_constrained=False)
                if N.all(g1 < 0): continue
            except DeadElement: continue
            (l2,g2) = p2.permuteElement(el2f, remove_constrained=False)
            if N.all(g2 < 0):
                if ignore_missing_2form: continue
                else: raise Exception(
                    'corresponding 2-form function is constrained')
            dofnos_1_g = g1[dofnos_1]
            fdofs = dofnos_1_g >= 0
            dofnos_2_g = g2[dofnos_2]
            if N.any(dofnos_2_g[fdofs] < 0) :
                if ignore_missing_2form: fdofs &= dofnos_2_g >= 0
                else: raise Exception(
                    'corresponding 2-form function is constrained')
            for d_2, d_1, w in izip(
                dofnos_2_g[fdofs], dofnos_1_g[fdofs], weights[fdofs]):
                assert(d_2 >= 0 and d_1 >= 0)
                mat[d_2][d_1] = w

        import array as A
        rows, cols, data = A.array('l'), A.array('l'), A.array('d')

        print 'Constructing row/column data'
        for i in xrange(len(mat)):
            row = mat[i]
            data.extend(row.values())
            cols_i = row.keys()
            del(row) ; mat[i] = None
            rows.extend([i]*len(cols_i))
            cols.extend(cols_i)
        del(mat)
        print 'Garbage Collecting'
        gc.collect()
        print 'Converting matrix formats'
        return Matrices.sparseType(sparse.coo_matrix(
            (data, (rows, cols)), dtype=self.dtype, shape=(no_2fdofs, no_1fdofs)))

class D_oneformTet(D_oneformBase):
    oneformType = 'webb99'
    twoformType = 'mmbotha06_solh'
    dtype = N.int8
    # This localmap hard codes the sorted local node numbering convention.
    curlLocalMaps = {(0.5, 0.5):N.array([[0, 0, 1], [0, 1, 1],
                                         [1, 0, -1], [1, 2, 1],
                                         [2, 1, -1], [2, 2, -1],
                                         [3, 0, 1], [3, 3, 1],
                                         [4, 1, 1], [4, 3, -1],
                                         [5, 2, 1], [5, 3, 1]], N.int8)}
    curlLocalMaps[(1.5, 1)] = N.vstack((curlLocalMaps[(0.5, 0.5)],
                                        N.array([[12,4,1], [13,5,1],
                                                 [14,6,1], [15,7,1],
                                                 [16,8,1], [17,9,1],
                                                 [18,10,1], [19,11,1]], N.int8)))
    curlLocalMaps[(2.5, 2)] = N.array([[0, 1, 1], [0, 0, 1],
                                       [1, 2, 1], [1, 0, -1],
                                       [2, 2, -1], [2, 1, -1],
                                       [3, 3, 1], [3, 0, 1],
                                       [4, 3, -1], [4, 1, 1],
                                       [5, 3, 1], [5, 2, 1],
                                       [18, 4, 1], [19, 5, 1],
                                       [20, 6, 1], [21, 7, 1],
                                       [22, 8, 1], [23, 9, 1],
                                       [24, 10, 1], [25, 11, 1],
                                       [30, 12, 1], [31, 13, 1],
                                       [32, 14, 1], [33, 15, 1],
                                       [34, 16, 1], [35, 17, 1],
                                       [36, 18, 1], [37, 19, 1],
                                       [38, 20, 1], [39, 21, 1],
                                       [40, 22, 1], [41, 23, 1],
                                       [42, 24, 1], [43, 25, 1],
                                       [44, 26, 1]], N.int8)
    curlLocalMaps[(3.5, 3)] = N.array([[0, 1, 1], [0, 0, 1],
                                       [1, 2, 1], [1, 0, -1],
                                       [2, 2, -1], [2, 1, -1],
                                       [3, 3, 1], [3, 0, 1],
                                       [4, 3, -1], [4, 1, 1],
                                       [5, 3, 1], [5, 2, 1],
                                       [24, 4, 1], [25, 5, 1],
                                       [26, 6, 1], [27, 7, 1],
                                       [28, 8, 1], [29, 9, 1],
                                       [30, 10, 1], [31, 11, 1],
                                       [36, 12, 1], [37, 13, 1],
                                       [38, 14, 1], [39, 15, 1],
                                       [40, 16, 1], [41, 17, 1],
                                       [42, 18, 1], [43, 19, 1],
                                       [44, 20, 1], [45, 21, 1],
                                       [46, 22, 1], [47, 23, 1],
                                       [56, 24, 1], [57, 25, 1],
                                       [58, 26, 1], [59, 27, 1],
                                       [60, 28, 1], [61, 29, 1],
                                       [62, 30, 1], [63, 31, 1],
                                       [64, 32, 1], [65, 33, 1],
                                       [66, 34, 1], [67, 35, 1],
                                       [68, 36, 1], [69, 37, 1],
                                       [70, 38, 1], [71, 39, 1],
                                       [72, 40, 1], [73, 41, 1],
                                       [74, 42, 1], [76, 43, 1],
                                       [77, 44, 1], [78, 45, 1],
                                       [79, 46, 1], [80, 47, 1],
                                       [81, 48, 1], [82, 49, 1],
                                       [83, 50, 1]], N.int8)


class D_oneformBrick(D_oneformBase):
    dtype = N.float64
    def _test_basis_types(self, disc1, disc2):
        from NewCode.Meshes import BrickMesh
        assert(isinstance(disc1.mesh, BrickMesh.Mesh))
        assert(isinstance(disc2.mesh, BrickMesh.Mesh))

    def _get_localMap(self, disc1, disc2):
        from NewCode.SystemMatrix import local_projection_matrix\
             ,local_self_projection_matrix, eps
        el_1 = disc1.D().elements[0]
        el_2 = disc2.elements[0]

        M_2_inv = linalg.inv(local_self_projection_matrix(el_2))
        P_21 = local_projection_matrix(el_2, el_1)
        C = N.dot(M_2_inv, P_21)
        C_abs = N.abs(C)
        C[C_abs/N.max(C_abs) < eps] = 0
        C = sparse.coo_matrix(C)
        return C.col, C.row, C.data.astype(self.dtype)

D_oneform_tet = D_oneformTet()
D_oneform_brick = D_oneformBrick()
