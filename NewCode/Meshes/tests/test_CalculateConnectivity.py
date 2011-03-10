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
        cc = self.calculate_connectivity
        cc.set_input_listmesh(self.input_listmesh)
        cc.setup_mesh()

    def test_setup_mesh(self):
        cc = self.calculate_connectivity
        assert_equal(cc.get_nodes(), self.desired_mesh.listmesh['Nodes'])
        assert_equal(cc.get_tet_nodes(), self.desired_mesh.listmesh['ElementNodes'])

    def test_calc_node_element_connectivity(self):
        cc = self.calculate_connectivity
        cc.calc_node_element_connectivity()
        assert_equal(cc.get_node_connect_2_element(),
                     self.desired_mesh.listmesh['NodeConnect2Element'])
        assert_equal(cc.get_node_connect_2_element_ptr(),
                     self.desired_mesh.listmesh['NodeConnect2ElementPtr'])

    def test_calc_edge_node_connectivity(self):
        cc = self.calculate_connectivity
        cc.calc_edge_node_connectivity()
        # Since the actual global edge numbering is arbitrary, we
        # check only check that the same set nof node pairs are
        # present
        desired_edge_nodes = set(
            tuple(en) for en in self.desired_mesh.listmesh['EdgeNodes'])
        e2n = cc.get_edge_connect_2_node()
        assert_equal(set(tuple(en) for en in e2n), desired_edge_nodes)


    def test_calc_element_edge_connectivity(self):
        cc = self.calculate_connectivity
        cc.calc_element_edge_connectivity()
        ref_edge_nodes = self.desired_mesh.listmesh['EdgeNodes']
        desired_element_edge_nodes = [
            ref_edge_nodes[i] for i in self.desired_mesh.listmesh['ElementEdges']]
        actual_element_connect_2_edge = cc.get_element_connect_2_edge()
        cc.calc_edge_node_connectivity()
        actual_edge_connect_2_node = cc.get_edge_connect_2_node()
        actual_element_connect_2_edge_nodes = array([
            actual_edge_connect_2_node[i] for i in actual_element_connect_2_edge])

        assert_equal(actual_element_connect_2_edge_nodes,
                     desired_element_edge_nodes)

    def test_calc_face_node_connectivity(self):
        cc = self.calculate_connectivity
        cc.calc_face_node_connectivity()
        # Since the actual global face numbering is arbitrary, we only
        # check that the same set nof node pairs are present
        desired_face_nodes = set(
            tuple(en) for en in self.desired_mesh.listmesh['FaceNodes'])
        f2n = cc.get_face_connect_2_node()
        assert_equal(set(tuple(fn) for fn in f2n), desired_face_nodes)

    def test_calc_element_face_connectivity(self):
        cc = self.calculate_connectivity
        cc.calc_element_face_connectivity()
        ref_face_nodes = self.desired_mesh.listmesh['FaceNodes']
        desired_element_face_nodes = [
            ref_face_nodes[i] for i in self.desired_mesh.listmesh['ElementFaces']]
        actual_element_connect_2_face = cc.get_element_connect_2_face()
        cc.calc_face_node_connectivity()
        actual_face_connect_2_node = cc.get_face_connect_2_node()
        actual_element_connect_2_face_nodes = array([
            actual_face_connect_2_node[i] for i in actual_element_connect_2_face])

        assert_equal(actual_element_connect_2_face_nodes,
                     desired_element_face_nodes)

    def test_calc_face_element_connectivity(self):
        cc = self.calculate_connectivity
        cc.calc_face_element_connectivity()
        ref_face_nodes = self.desired_mesh.listmesh['FaceNodes']
        ref_face_elements = self.desired_mesh.listmesh['FaceConnect2Elem']
        desired_facenode_elements = dict((tuple(facenodes), els) for facenodes, els in zip(
            ref_face_nodes, ref_face_elements))
        actual_face_elements = cc.get_face_connect_2_element()
        cc.calc_face_node_connectivity()
        actual_face_nodes = cc.get_face_connect_2_node()
        actual_facenode_elements = dict((tuple(facenodes), els) for facenodes, els in zip(
            actual_face_nodes, actual_face_elements))
        assert_equal(actual_facenode_elements, desired_facenode_elements)
