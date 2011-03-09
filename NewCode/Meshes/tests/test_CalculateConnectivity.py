from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, \
     assert_equal, assert_array_almost_equal
from numpy import array, float64, int32
import sys, os
#
# Local Imports
#
from NewCode.Meshes import CalculateConnectivity
import NewCode.tests.TestMeshes as TestMeshes

InscribedTetMesh = TestMeshes.InscribedTetMesh


class test_CalculateConnectivity(TestCase):
    DesiredMesh = InscribedTetMesh
    def setUp(self):
        self.desired_mesh = self.DesiredMesh()
        self.desired_mesh.setUp()
        lm = self.desired_mesh.listmesh
        self.input_listmesh = dict(Nodes=lm['Nodes'].copy(),
                                   ElementNodes=lm['ElementNodes'].copy())
        self.calculate_connectivity = CalculateConnectivity.CalculateConnectivity()
        
    
    def test_setup_mesh(self):
        cc = self.calculate_connectivity
        cc.set_listmesh(self.input_listmesh)
        cc.setup_mesh()
        assert_equal(cc.get_nodes(), self.desired_mesh.listmesh['Nodes'])
        assert_equal(cc.get_tet_nodes(), self.desired_mesh.listmesh['ElementNodes'])
