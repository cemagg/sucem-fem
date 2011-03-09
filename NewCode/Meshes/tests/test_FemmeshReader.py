from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, \
     assert_equal, assert_array_almost_equal
from numpy import array, float64, int32
import sys, os
#
# Local Imports
#
from NewCode.Meshes import FemmeshReader
import NewCode.tests.TestMeshes as TestMeshes

InscribedTetMesh = TestMeshes.InscribedTetMesh


class test_FemmeshReader(TestCase):
    DesiredMesh = InscribedTetMesh
    femmesh_filename = 'NewCode/tests/testdata/test_FemmeshReader/inscribed_tet.femmesh'
    def setUp(self):
        self.desired_mesh = self.DesiredMesh()
        self.desired_mesh.setUp()
        self.femmesh_reader = FemmeshReader.FemmeshReader(self.femmesh_filename)
        pass
    
    def test_read_meshfile(self):
        self.femmesh_reader.read_meshfile()
        assert_almost_equal(self.femmesh_reader.nodes,
                            self.desired_mesh.listmesh['Nodes'],
                            decimal=14)
        assert_equal(self.femmesh_reader.tet_nodes,
                     self.desired_mesh.listmesh['ElementNodes'])
        
    def test_parse_nodes(self):
        parse_str = ['           8\n',
                     '           1   -0.5                  0.5                      -0.5\n',
                     '           2   0.5                   0.5                      0.5\n',
                     '           3   0.5                   -0.5                     -0.5\n',
                     '           4   -0.5                  -0.5                     0.5\n',
                     '           5   0.1666666666666667     0.1666666666666667      -0.1666666666666667\n',
                     '           6   -0.1666666666666667    0.1666666666666667      0.1666666666666667\n',
                     '           7   -0.1666666666666667    -0.1666666666666667     -0.1666666666666667\n',
                     '           8   0.1666666666666667     -0.1666666666666667     0.1666666666666667\n',
                     'ENDBLOCK\n']
        self.femmesh_reader.parse_nodes(iter(parse_str))
        assert_almost_equal(self.femmesh_reader.nodes,
                            self.desired_mesh.listmesh['Nodes'],
                            decimal=14)


    def test_parse_tets(self):
        parse_str = ['          11\n',
        '           1           1           1           6           4           7\n',
        '           2           1           1           3           5           7\n',
        '           3           1           4           3           7           8\n',
        '           4           1           2           1           5           6\n',
        '           5           1           3           2           5           8\n',
        '           6           1           2           4           6           8\n',
        '           7           1           8           4           6           7\n',
        '           8           1           1           5           6           7\n',
        '           9           1           3           5           7           8\n',
        '          10           1           5           2           6           8\n',
        '          11           1           5           6           7           8\n',
        'ENDBLOCK\n']
        self.femmesh_reader.parse_tets(iter(parse_str))
        assert_equal(self.femmesh_reader.tet_nodes,
                     self.desired_mesh.listmesh['ElementNodes'])

class test_Femmesh2ListMesh(TestCase):
    DesiredMesh = InscribedTetMesh
    femmesh_filename = 'NewCode/tests/testdata/test_FemmeshReader/inscribed_tet.femmesh'
    def setUp(self):
        self.desired_mesh = self.DesiredMesh()
        self.desired_mesh.setUp()
        self.femmesh_reader = FemmeshReader.FemmeshReader(self.femmesh_filename)
        self.femmesh_reader.read_meshfile()
        self.femmesh2listmesh = FemmeshReader.Femmesh2ListMesh(self.femmesh_reader)
        
    def test_listmesh_filename(self):
        testlistmesh = self.femmesh2listmesh.get_listmesh()
        assert_equal(testlistmesh['FemmeshFilename'], 'inscribed_tet.femmesh')

    def test_listmesh_nodes(self):
        testlistmesh = self.femmesh2listmesh.get_listmesh()
        assert_almost_equal(testlistmesh['Nodes'], self.desired_mesh.listmesh['Nodes'],
                            decimal=14)

    def test_listmesh_elementnodes(self):
        testlistmesh = self.femmesh2listmesh.get_listmesh()
        assert_equal(testlistmesh['ElementNodes'], self.desired_mesh.listmesh['ElementNodes'])
