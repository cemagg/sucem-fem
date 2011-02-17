from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.DifferentialForm import Discretiser, DiscretiserEntities
from NewCode.DifferentialForm import DirechletDiscretiser
from NewCode.tests.TestMeshes import FlatTet, TwoTets, InscribedTetMesh
from NewCode import Mesh

class test_DirechletDiscretiser(NumpyTestCase):
    @staticmethod
    def freeFun(ent):
        ns = [1,2,3,5]
        return N.all([x+1 in ns for x in ent.nodes])
    
    def setUp(self):
        self.mesh = Mesh.Mesh(InscribedTetMesh.listmesh)
        self.disc = Discretiser.setup_PformDiscretiser(self.mesh, 1)
        self.geomEntities = DiscretiserEntities.make_geomEntities(
            self.mesh, self.disc.basisSet, self.freeFun)

    def test_init(self):
        inst = DirechletDiscretiser.DirechletDiscretiser(
            self.disc, self.geomEntities, Discretiser.Permuter)


    def test_matrix(self):
        inst = DirechletDiscretiser.DirechletDiscretiser(
            self.disc, self.geomEntities, Discretiser.Permuter)
        disc_mass = self.disc.matrix.mass().todense()
        direch_mat = inst.matrix.projectionOnto(self.disc).todense()
        assert_equal(direch_mat.shape,
                     (disc_mass.shape[0], self.geomEntities['edge'].noFree))
        assert_equal(direch_mat, disc_mass[:, [e.index for e in self.geomEntities['edge']
                                               if e.isFree]])

    def test_matrix_D(self):
        inst = DirechletDiscretiser.DirechletDiscretiser(
            self.disc, self.geomEntities, Discretiser.Permuter)
        disc_stiff = self.disc.matrix.stiffness().todense()
        direch_mat = inst.D().matrix.projectionOnto(self.disc.D()).todense()
        assert_equal(direch_mat.shape,
                     (disc_stiff.shape[0], self.geomEntities['edge'].noFree))
        assert_equal(direch_mat, disc_stiff[:, [e.index for e in self.geomEntities['edge']
                                               if e.isFree]])
        
