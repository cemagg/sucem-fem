from __future__ import division
from numpy.testing import TestCase, assert_array_equal, \
     assert_array_almost_equal, assert_almost_equal, assert_equal
import numpy as N
from scipy import sparse 

from NewCode.tests.TestMeshes import FlatTet, TwoTets
from NewCode.tests.PyramMeshes import SixPyram
from NewCode.tests import xfail
from NewCode.Utilities import Struct
from NewCode import SystemMatrix, Mesh, Integration
from NewCode.DifferentialForm import Discretiser, DiscretiserEntities
from NewCode.DifferentialForm import BrickDiscretiser, BrickDiscretiserEntities
from NewCode.DifferentialForm import PyramDiscretiser, PyramDiscretiserEntities
from NewCode.Meshes import BrickMesh, PyramMesh
import SystemMatrixValues, BrickSystemMatrixValues, PyramSystemMatrixValues

class test_matrices(TestCase):
    def setUp(self):
        self.mesh = Mesh.Mesh(TwoTets.listmesh)

    def _setup_oneform_twotet(self, order=1, mixed=True):
        from NewCode.DifferentialForm.BasisFunction import Oneform
        basisSet = Oneform.basis_set(order=order, mixed=mixed)
        geomEntities = DiscretiserEntities.make_geomEntities(
            self.mesh, basisSet, lambda x: True)
        return  Discretiser.PformDiscretiser(
            1, self.mesh, geomEntities, Discretiser.Permuter, basisSet)

    def _setup_twoform_twotet(self, order=1, mixed=True):
        from NewCode.DifferentialForm.BasisFunction import Twoform
        basisSet = Twoform.basis_set(order=order, mixed=mixed)
        geomEntities = DiscretiserEntities.make_geomEntities(
            self.mesh, basisSet, lambda x: True)
        return  Discretiser.PformDiscretiser(
            2, self.mesh, geomEntities, Discretiser.Permuter, basisSet)
    
class test_self_projection_matrix_oneform(test_matrices):
    def test(self):
        mat = SystemMatrix.self_projection_matrix(Struct(
            totalDOFs=5, elements=[], identicalElements=False,
            diagonalMetric=False))
        assert_equal(mat.shape, (5,5))
        assert_equal(mat.dtype, N.dtype(N.float64))
        
    def test_value_oneform(self):
        disc = self._setup_oneform_twotet()
        mat = SystemMatrix.self_projection_matrix(disc)
        desired_mat = SystemMatrixValues.global_onef_1_mixed_mass
        actual_mat = mat.toarray()
        assert_array_almost_equal(actual_mat, desired_mat, decimal=15)

    def test_value_oneform_localmats(self):
        disc = self._setup_oneform_twotet()
        el0 = disc.elements[0]
        lm0 = SystemMatrix.local_self_projection_matrix(el0)
        el1 = disc.elements[1]
        lm1 = SystemMatrix.local_self_projection_matrix(el1)
        lm0_des, lm1_des = SystemMatrixValues.local_onef_1_mixed_mass
        assert_array_almost_equal(lm0, lm0_des, decimal=15)
        assert_array_almost_equal(lm1, lm1_des, decimal=15)

    def test_value_oneform_LTQN(self):
        disc = self._setup_oneform_twotet(2)
        mat = SystemMatrix.self_projection_matrix(disc)
        desired_mat = SystemMatrixValues.global_onef_2_mixed_mass
        actual_mat = mat.toarray()
        assert_array_almost_equal(actual_mat, desired_mat, decimal=15)
        
    def test_value_oneform_localmats_LTQN(self):
        disc = self._setup_oneform_twotet(order=2)
        el0 = disc.elements[0]
        lm0 = SystemMatrix.local_self_projection_matrix(el0)
        el1 = disc.elements[1]
        lm1 = SystemMatrix.local_self_projection_matrix(el1)
        lm0_des, lm1_des = SystemMatrixValues.local_onef_2_mixed_mass
        assert_array_almost_equal(lm0, lm0_des, decimal=15)
        assert_array_almost_equal(lm1, lm1_des, decimal=15)

    def test_value_oneform_QTQN(self):
        disc = self._setup_oneform_twotet(2, mixed=False)
        mat = SystemMatrix.self_projection_matrix(disc)
        desired_mat = SystemMatrixValues.global_onef_2_mass
        actual_mat = mat.toarray()
        assert_array_almost_equal(actual_mat, desired_mat, decimal=15)
        
    def test_value_oneform_localmats_QTQN(self):
        disc = self._setup_oneform_twotet(order=2, mixed=False)
        el0 = disc.elements[0]
        lm0 = SystemMatrix.local_self_projection_matrix(el0)
        el1 = disc.elements[1]
        lm1 = SystemMatrix.local_self_projection_matrix(el1)
        lm0_des, lm1_des = SystemMatrixValues.local_onef_2_mass
         
        assert_array_almost_equal(lm0, lm0_des, decimal=15)
    
    def test_value_oneform_4th(self):
        disc = self._setup_oneform_twotet(4, mixed=False)
        mat = SystemMatrix.self_projection_matrix(disc)
        desired_mat = SystemMatrixValues.global_onef_4_mass
        actual_mat = mat.toarray()
        assert_array_almost_equal(actual_mat, desired_mat, decimal=15)
        
    def test_value_oneform_localmats_4th(self):
        disc = self._setup_oneform_twotet(order=4, mixed=False)
        el0 = disc.elements[0]
        lm0 = SystemMatrix.local_self_projection_matrix(el0)
        el1 = disc.elements[1]
        lm1 = SystemMatrix.local_self_projection_matrix(el1)
        lm0_des, lm1_des = SystemMatrixValues.local_onef_4_mass
        
        assert_array_almost_equal(lm0, lm0_des, decimal=15)
        assert_array_almost_equal(lm1, lm1_des, decimal=15)

    def test_value_oneform_D(self):
        from NewCode.DifferentialForm.BasisFunction import Oneform
        mesh = Mesh.Mesh(TwoTets.listmesh)
        basisSet = Oneform.basis_set(1)
        geomEntities = DiscretiserEntities.make_geomEntities(mesh, basisSet, lambda x: True)
        disc = Discretiser.PformDiscretiser(
            1, mesh, geomEntities, Discretiser.Permuter, basisSet)
        mat = SystemMatrix.self_projection_matrix(disc.D()).todense()
        desired_mat = SystemMatrixValues.global_onef_1_mixed_stiffness
        assert_array_almost_equal(mat, desired_mat, decimal=14)
        
    def test_value_oneform_localmats_D(self):
        disc = self._setup_oneform_twotet().D()
        
        el0 = disc.elements[0]
        lm0 = SystemMatrix.local_self_projection_matrix(el0)
        el1 = disc.elements[1]
        lm1 = SystemMatrix.local_self_projection_matrix(el1)
        lm0_des, lm1_des = SystemMatrixValues.local_onef_1_mixed_stiffness
        assert_array_almost_equal(lm0, lm0_des, decimal=14)
        assert_array_almost_equal(lm1, lm1_des, decimal=14)
        
    def test_value_oneform_LTQN_D(self):
        disc = self._setup_oneform_twotet(2).D()
        mat = SystemMatrix.self_projection_matrix(disc)
        desired_mat = SystemMatrixValues.global_onef_2_mixed_stiffness
        actual_mat = mat.toarray()
        assert_array_almost_equal(actual_mat, desired_mat, decimal=14)

    def test_value_oneform_localmats_LTQN_D(self):
        disc = self._setup_oneform_twotet(order=2).D()
        el0 = disc.elements[0]
        lm0 = SystemMatrix.local_self_projection_matrix(el0)
        el1 = disc.elements[1]
        lm1 = SystemMatrix.local_self_projection_matrix(el1)
        lm0_des, lm1_des = SystemMatrixValues.local_onef_2_mixed_stiffness
        assert_array_almost_equal(lm0, lm0_des, decimal=14)
        assert_array_almost_equal(lm1, lm1_des, decimal=14)

    def test_value_oneform_QTQN_D(self):
        disc = self._setup_oneform_twotet(2, mixed=False).D()
        mat = SystemMatrix.self_projection_matrix(disc)
        desired_mat = SystemMatrixValues.global_onef_2_stiffness
        actual_mat = mat.toarray()
        assert_array_almost_equal(actual_mat, desired_mat, decimal=14)

    def test_value_oneform_localmats_QTQN_D(self):
        disc = self._setup_oneform_twotet(order=2, mixed=False).D()
        el0 = disc.elements[0]
        lm0 = SystemMatrix.local_self_projection_matrix(el0)
        el1 = disc.elements[1]
        lm1 = SystemMatrix.local_self_projection_matrix(el1)
        lm0_des,lm1_des = SystemMatrixValues.local_onef_2_stiffness
        assert_array_almost_equal(lm0, lm0_des, decimal=14)

    def test_value_oneform_4th_D(self):
        disc = self._setup_oneform_twotet(order=4, mixed=False).D()
        mat = SystemMatrix.self_projection_matrix(disc)
        desired_mat = SystemMatrixValues.global_onef_4_stiffness
        actual_mat = mat.toarray()
        assert_array_almost_equal(actual_mat, desired_mat, decimal=14)

    def test_value_oneform_localmats_4th_D(self):
        disc = self._setup_oneform_twotet(order=4, mixed=False).D()
        el0 = disc.elements[0]
        lm0 = SystemMatrix.local_self_projection_matrix(el0)
        el1 = disc.elements[1]
        lm1 = SystemMatrix.local_self_projection_matrix(el1)
        lm0_des, lm1_des = SystemMatrixValues.local_onef_4_stiffness
        assert_array_almost_equal(lm0, lm0_des, decimal=14)
        assert_array_almost_equal(lm1, lm1_des, decimal=14)

class test_self_projection_matrix_twoform(test_matrices):
    def test_value_twoform_1_mixed(self):
        disc = self._setup_twoform_twotet()
        desired_mat = SystemMatrixValues.global_twof_1_mixed_mass
        actual_mat = SystemMatrix.self_projection_matrix(disc)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

    def test_value_twoform_1(self):
        disc = self._setup_twoform_twotet(mixed=False)
        desired_mat = SystemMatrixValues.global_twof_1_mass
        actual_mat = SystemMatrix.self_projection_matrix(disc)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)
        desired_mat_D = SystemMatrixValues.global_twof_1_stiffness
        actual_mat_D = SystemMatrix.self_projection_matrix(disc.D())
        assert_almost_equal(actual_mat_D.toarray(), desired_mat_D, decimal=13)

    def test_value_twoform_1_mixed_D(self):
        pass


    
class test_projection_matrix(test_matrices):
    def test_value_oneform_twoform(self):
        desired_mat = SystemMatrixValues.global_onef_1_mixed_twof_1_mixed_proj
        from NewCode.DifferentialForm.BasisFunction import Oneform, Twoform
        disc_oneform = self._setup_oneform_twotet()
        disc_twoform = self._setup_twoform_twotet()
        assert_almost_equal(SystemMatrix.projection_matrix(
            disc_oneform, disc_twoform).toarray(),
                            desired_mat, decimal=14)

class test_boundary_matrix(test_matrices):
    def test_value_oneform_CTLN_localmats(self):
        des_lm_0, des_lm_1 = SystemMatrixValues.local_onef_1_mixed_boundary_surface_integral
        disc = self._setup_oneform_twotet(order=1, mixed=True)
        mesh = Mesh.Mesh(TwoTets.listmesh)
        disc.elements[:].setFaceIntegrationRule(Integration.TriIntegrator(2))
        el0 = disc.elements[0]
        lm0 = SystemMatrix.local_boundary_matrix(el0, disc.mesh.faces)
        assert_array_almost_equal(lm0, des_lm_0, decimal=16)
        el1 = disc.elements[1]
        lm1 = SystemMatrix.local_boundary_matrix(el1, disc.mesh.faces)
        assert_array_almost_equal(lm1, des_lm_1, decimal=16)

    def test_value_oneform_mixed2_localmats(self):
        des_lm_0, des_lm_1 = SystemMatrixValues.local_onef_2_mixed_boundary_surface_integral
        disc = self._setup_oneform_twotet(order=2, mixed=True)
        mesh = Mesh.Mesh(TwoTets.listmesh)
        disc.elements[:].setFaceIntegrationRule(Integration.TriIntegrator(4))
        el0 = disc.elements[0]
        lm0 = SystemMatrix.local_boundary_matrix(el0, disc.mesh.faces)
        assert_array_almost_equal(lm0, des_lm_0, decimal=15)
        el1 = disc.elements[1]
        lm1 = SystemMatrix.local_boundary_matrix(el1, disc.mesh.faces)
        assert_array_almost_equal(lm1, des_lm_1, decimal=15)
        

class test_insert_global(TestCase):
    def test_one_index(self):
        # Permutation for "element" 1
        perm_el1 = (N.arange(2), N.array([2,1], int))
        # Permutation for "element" 2
        perm_el2 = (N.arange(2), N.array([1,0], int))
        # local matrix values
        val_el1 = N.array([[10,100],
                           [1000, 10000]], int)
        val_el2 = N.array([[15, 150],
                           [1500, 15000]], int)
        glob_mat = sparse.lil_matrix((3,3), dtype=int)
        SystemMatrix.insert_global(glob_mat, val_el1, perm_el1)
        SystemMatrix.insert_global(glob_mat, val_el2, perm_el2)
        desired = N.array([[15000, 1500,0],
                           [150, 10015,1000],
                           [0,100,10]])
        assert_equal(glob_mat.todense(), desired)

    def test_two_indices(self):
        pass

class _test_BrickMatrices(TestCase):
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.mesh = BrickMesh.Mesh(self.testMesh.listmesh)

    def _setup_oneform(self, order=1, mixed=True, btype=None):
        from NewCode.DifferentialForm.BasisFunction import BrickOneform
        if btype:
            basisSet = BrickOneform.basis_set(order=order, mixed=mixed, btype=btype)
        else:
            basisSet = BrickOneform.basis_set(order=order, mixed=mixed)
        geomEntities = BrickDiscretiserEntities.make_geomEntities(
            self.mesh, basisSet, lambda x: True)
        return  BrickDiscretiser.PformDiscretiser(
            1, self.mesh, geomEntities, Discretiser.Permuter, basisSet)

    def _setup_twoform(self, order=1, mixed=True, btype=None):
        from NewCode.DifferentialForm.BasisFunction import BrickTwoform
        if btype:
            basisSet = BrickTwoform.basis_set(order=order, mixed=mixed, btype=btype)
        else:
            basisSet = BrickTwoform.basis_set(order=order, mixed=mixed)
        geomEntities = BrickDiscretiserEntities.make_geomEntities(
            self.mesh, basisSet, lambda x: True)
        return  BrickDiscretiser.PformDiscretiser(
            2, self.mesh, geomEntities, Discretiser.Permuter, basisSet)

class test_BrickMatrices_OneBrick(_test_BrickMatrices):
    from NewCode.tests.BrickMeshes import OneBrick as TestMesh

    def test_value_oneform_mixed1_localmat(self):
        disc = self._setup_oneform()
        el0 = disc.elements[0]
        lm0 = SystemMatrix.local_self_projection_matrix(el0)
        assert_almost_equal(lm0, BrickSystemMatrixValues.local_onef_1_mixed_mass,
                            decimal=14)

        el0 = disc.D().elements[0]
        lm0_D = SystemMatrix.local_self_projection_matrix(el0)
        assert_almost_equal(lm0_D, BrickSystemMatrixValues.local_onef_1_mixed_stiffness,
                            decimal=15)


class test_BrickMatrices_TwoBrick(_test_BrickMatrices):
    from NewCode.tests.BrickMeshes import TwoBricks as TestMesh

    def test_value_oneform_1_mixed(self):
        disc = self._setup_oneform()
        desired_mat = BrickSystemMatrixValues.global_onef_1_mixed_mass
        actual_mat = SystemMatrix.self_projection_matrix(disc)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=15)

        desired_mat = BrickSystemMatrixValues.global_onef_1_mixed_stiffness
        actual_mat = SystemMatrix.self_projection_matrix(disc.D())
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

    def test_value_2_mixed_rieben04(self):
        disc = self._setup_oneform(order=2)
        disc2 = self._setup_twoform(order=2)
        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_2_mixed_rieben04_mass'].todense()
        actual_mat = SystemMatrix.self_projection_matrix(disc)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_2_mixed_rieben04_stiffness'].todense()
        actual_mat = SystemMatrix.self_projection_matrix(disc.D())
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_twoform_2_mixed_rieben04_P'].todense()
        actual_mat = SystemMatrix.projection_matrix(disc.D(), disc2)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)
        
    def test_value_4th_mixed_rieben04(self):
        disc = self._setup_oneform(order=4)
        disc2 = self._setup_twoform(order=4)
        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_4_mixed_rieben04_mass'].todense()
        actual_mat = SystemMatrix.self_projection_matrix(disc)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_4_mixed_rieben04_stiffness'].todense()
        actual_mat = SystemMatrix.self_projection_matrix(disc.D())
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_twoform_4_mixed_rieben04_P'].todense()
        actual_mat = SystemMatrix.projection_matrix(disc.D(), disc2)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)
        
    def test_value_2_mixed_cohen98(self):
        disc = self._setup_oneform(order=2, btype='cohen98')
        disc2 = self._setup_twoform(order=2, btype='cohen98')
        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_2_mixed_cohen98_mass'].todense()
        actual_mat = SystemMatrix.self_projection_matrix(disc)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_2_mixed_cohen98_stiffness'].todense()
        actual_mat = SystemMatrix.self_projection_matrix(disc.D())
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_twoform_2_mixed_cohen98_P'].todense()
        actual_mat = SystemMatrix.projection_matrix(disc.D(), disc2)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

    def test_value_4th_mixed_cohen98(self):
        disc = self._setup_oneform(order=4, btype='cohen98')
        disc2 = self._setup_twoform(order=4, btype='cohen98')
        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_4_mixed_cohen98_mass'].todense()
        actual_mat = SystemMatrix.self_projection_matrix(disc)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_4_mixed_cohen98_stiffness'].todense()
        actual_mat = SystemMatrix.self_projection_matrix(disc.D())
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=13)

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_twoform_4_mixed_cohen98_P'].todense()
        actual_mat = SystemMatrix.projection_matrix(disc.D(), disc2)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)
        
    def test_value_4th_mixed_cohen98_lumped(self):
        order = 4
        disc = self._setup_oneform(order=order, btype='cohen98')
        disc2 = self._setup_twoform(order=order, btype='cohen98')
        disc.diagonalMetric = True
        disc2.diagonalMetric = True
        disc.diagonalise()
        disc2.diagonalise()

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_4_mixed_cohen98_lumped_mass'].todense()
        actual_mat = SystemMatrix.self_projection_matrix(disc)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_4_mixed_cohen98_lumped_stiffness'].todense()
        actual_mat = SystemMatrix.self_projection_matrix(disc.D())
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_twoform_4_mixed_cohen98_lumped_P'].todense()
        actual_mat = SystemMatrix.projection_matrix(disc.D(), disc2)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)
        
    def test_value_2_mixed_cohen98_lumped(self):
        order = 2
        disc = self._setup_oneform(order=order, btype='cohen98')
        disc2 = self._setup_twoform(order=order, btype='cohen98')
        disc.diagonalMetric = True
        disc2.diagonalMetric = True
        disc.diagonalise()
        disc2.diagonalise()

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_2_mixed_cohen98_lumped_mass'].todense()
        actual_mat = SystemMatrix.self_projection_matrix(disc)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_2_mixed_cohen98_lumped_stiffness'].todense()
        actual_mat = SystemMatrix.self_projection_matrix(disc.D())
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)

        desired_mat = BrickSystemMatrixValues.more_brick_mats[
            'global_oneform_twoform_2_mixed_cohen98_lumped_P'].todense()
        actual_mat = SystemMatrix.projection_matrix(disc.D(), disc2)
        assert_almost_equal(actual_mat.toarray(), desired_mat, decimal=14)
        
class _test_PyramMatrices(TestCase):
    TestMesh = SixPyram
    
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.mesh = PyramMesh.Mesh(self.testMesh.listmesh)

    def _setup_oneform(self, order=1, mixed=True, btype=None):
        from NewCode.DifferentialForm.BasisFunction import PyramOneform
        if btype:
            basisSet = PyramOneform.basis_set(order=order, mixed=mixed, btype=btype)
        else:
            basisSet = PyramOneform.basis_set(order=order, mixed=mixed)
        geomEntities = PyramDiscretiserEntities.make_geomEntities(
            self.mesh, basisSet, lambda x: True)
        return  PyramDiscretiser.PformDiscretiser(
            1, self.mesh, geomEntities, Discretiser.Permuter, basisSet)

class test_PyramOneform(_test_PyramMatrices):
    def test_local_mixed_1_mass(self):
        disc = self._setup_oneform(btype='graglia99')
        desired_mats = PyramSystemMatrixValues.local_onef_1_mixed_mass
        assert_almost_equal([
            SystemMatrix.local_self_projection_matrix(el)
            for el in disc.elements],
                            N.array(desired_mats), decimal=15)

    def test_local_mixed_1_stiffness(self):
        disc = self._setup_oneform(btype='graglia99').D()
        desired_mats = PyramSystemMatrixValues.local_onef_1_mixed_stiffness
        assert_almost_equal([
            SystemMatrix.local_self_projection_matrix(el)
            for el in disc.elements],
                            N.array(desired_mats), decimal=15)
