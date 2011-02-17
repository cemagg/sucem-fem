from __future__ import division
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal
import numpy as N
from scipy import sparse 

from NewCode.tests.TestMeshes import FlatTet, TwoTets
from NewCode.tests import xfail
from NewCode.Utilities import Struct
from NewCode import SystemMatrix, Mesh, SubDimMesh
from NewCode.DifferentialForm import Discretiser, DiscretiserEntities
from NewCode.DifferentialForm import SubDimDiscretiser, SubDimDiscretiserEntities

class test_matrices(NumpyTestCase):
    faceSelector = staticmethod(lambda face: 1 in face.connect2elem)
    freefun = staticmethod(lambda x: True)
    def setUp(self):
        self.superMesh = Mesh.Mesh(TwoTets.listmesh)
        self.mesh = SubDimMesh.SubSurface(self.superMesh, self.faceSelector)

    def _setup_oneform_superdisc(self, order=1, mixed=True):
        from NewCode.DifferentialForm.BasisFunction import Oneform
        basisSet = Oneform.basis_set(order=order, mixed=mixed)
        geomEntities = DiscretiserEntities.make_geomEntities(
            self.superMesh, basisSet, lambda x: True)
        return  Discretiser.PformDiscretiser(
            1, self.superMesh, geomEntities, Discretiser.Permuter, basisSet)

    def _setup_oneform_disc(self, mesh, superMesh, superDisc, freefun=None):
        if freefun is None: freefun = self.freefun
        geomEntities = {'edge': SubDimDiscretiserEntities.SubDimDiscretiserEntityList(
            SubDimDiscretiserEntities.Edge(mesh=mesh, freefun=freefun,
                                           attrs=mesh.edges.list_repr())),
                        }
        return SubDimDiscretiser.PformSubDimDiscretiser(
           1, mesh, geomEntities, Discretiser.Permuter, superDisc)

#### NOTE!!! No values are tested yet, just looking for "syntax errors" #####

    def test_value_mass(self):
        superDisc = self._setup_oneform_superdisc()
        disc = self._setup_oneform_disc(self.mesh, self.superMesh, superDisc)
        mat = SystemMatrix.self_projection_matrix(disc)
        print mat.todense()

    def test_value_stiffness(self):
        superDisc = self._setup_oneform_superdisc()
        disc = self._setup_oneform_disc(self.mesh, self.superMesh, superDisc)
        mat = SystemMatrix.self_projection_matrix(disc.D())
        print mat.todense()
