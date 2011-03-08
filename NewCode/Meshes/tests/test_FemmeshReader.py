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
        #self.testlistmesh = self.femmesh_reader.get_listmesh()
        pass
    
    # def test_femmesh_filename(self):
    #     assert_equal(self.testlistmesh['FemmeshFilename'], 'inscribed_tet.femmesh')

    # def test_femmesh_nodes(self):
    #     assert_almost_equal(self.testlistmesh['Nodes'],
    #                         self.desired_mesh.listmesh['Nodes'])

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


