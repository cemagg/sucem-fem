from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal
from NewCode.Utilities import Struct
from NewCode import DiscretisedSystem as DS
from NewCode import Mesh, DifferentialForm
from NewCode.DifferentialForm import Discretiser
from NewCode.tests.TestMeshes import TwoTets, InscribedTetMesh

class test_VectorWaveEigen(NumpyTestCase):
    def setUp(self):
        self.mesh = Mesh.Mesh(TwoTets.listmesh)
        self.eigsys=DS.VectorWaveEigen(self.mesh,1, BC=lambda *x: True)

    def test_massMatrix(self):
        desired_mat = N.array([[7/120,1/80,1/80,-1/80,-1/80,0,0,0,0],
                               [1/80,11/120,1/48,1/80,0,-1/48,0,0,0],
                               [1/80,1/48,11/120,0,1/80,1/48,0,0,0],
                               [-1/80,1/80,0,7/120,1/80,-1/80,0,0,0],
                               [-1/80,0,1/80,1/80,7/120,1/80,0,0,0],
                               [0,-1/48,1/48,-1/80,1/80,11/120,0,0,0],
                               [0,0,0,0,0,0,1/12,1/24,1/24],
                               [0,0,0,0,0,0,1/24,1/12,1/24],
                               [0,0,0,0,0,0,1/24,1/24,1/12]], N.float64)
        calc = self.eigsys.massMatrix()
        assert_almost_equal(desired_mat,calc.todense(),
                            decimal=14)
        self.assert_(calc is self.eigsys.massMatrix()) # test caching

    def test_stiffnessMatrix(self):
        desired_mat = N.array([[2/3,-1/3,-1/3,1/3,1/3,0,0,0,0],
                               [-1/3,4/3,-1/3,-1/3,0,1/3,-2/3,2/3,0],
                               [-1/3,-1/3,4/3,0,-1/3,-1/3,-2/3,0,2/3],
                               [1/3,-1/3,0,2/3,-1/3,1/3,0,0,0],
                               [1/3,0,-1/3,-1/3,2/3,-1/3,0,0,0],
                               [0,1/3,-1/3,1/3,-1/3,4/3,0,-2/3,2/3],
                               [0,-2/3,-2/3,0,0,0,4/3,-2/3,-2/3],
                               [0,2/3,0,0,0,-2/3,-2/3,4/3,-2/3],
                               [0,0,2/3,0,0,2/3,-2/3,-2/3,4/3]], N.float64)
        calc = self.eigsys.stiffnessMatrix()
        assert_almost_equal(desired_mat, calc.todense(),
                            decimal=14)
        self.assert_(calc is self.eigsys.stiffnessMatrix()) # test caching
        

class test_AudioEigen(NumpyTestCase):
    def setUp(self):
        self.mesh = Mesh.Mesh(TwoTets.listmesh)
        self.eigsys=DS.AudioEigen(self.mesh,1, BC=lambda *x: True)

    def test_massMatrix(self):
        desired_mat = N.array([[3/10,1/30,-1/30,1/30,0,0,0],
                               [1/30,3/10,1/30,-1/30,0,0,0],
                               [-1/30,1/30,1/2,1/30,-1/30,1/30,-1/30],
                               [1/30,-1/30,1/30,3/10,0,0,0],
                               [0,0,-1/30,0,8/15,2/15,-2/15],
                               [0,0,1/30,0,2/15,8/15,2/15],
                               [0,0,-1/30,0,-2/15,2/15,8/15]], N.float64)
        calc = self.eigsys.massMatrix()
        assert_almost_equal(desired_mat,calc.todense(),
                            decimal=14)
        self.assert_(calc is self.eigsys.massMatrix()) # test caching
    def test_stiffnessMatrix(self):
        desired_mat = N.array([[3,-3,3,-3,0,0,0],
                               [-3,3,-3,3,0,0,0],
                               [3,-3,9,-3,-6,6,-6],
                               [-3,3,-3,3,0,0,0],
                               [0,0,-6,0,6,-6,6],
                               [0,0,6,0,-6,6,-6],
                               [0,0,-6,0,6,-6,6]], N.float64)
        calc = self.eigsys.stiffnessMatrix()
        assert_almost_equal(desired_mat, calc.todense(),
                            decimal=13)
        self.assert_(calc is self.eigsys.stiffnessMatrix()) # test caching

class test_CoupledFirstOrderSystem(NumpyTestCase):
    def setUp(self):
        self.mesh = mesh = Mesh.Mesh(InscribedTetMesh.listmesh)
        cb = DifferentialForm.constrained_on_boundary
        free = DifferentialForm.allfree        
        self.BCs = BCs = Struct(E=cb, H=free, D=free, B=cb)
        self.inst = DS.CoupledFirstOrderSystem(mesh, BCs)
        self.inst.setTimestep(0.5)

    def test_init(self):
        inst = self.inst ; BCs = self.BCs ; mesh = self.mesh
        assert_equal(inst.dt, 0.5)
        self.assert_(isinstance(inst.discs.E, Discretiser.PformDiscretiser))
        self.assert_(isinstance(inst.discs.D, Discretiser.PformDiscretiser))
        self.assert_(isinstance(inst.discs.H, Discretiser.PformDiscretiser))
        self.assert_(isinstance(inst.discs.B, Discretiser.PformDiscretiser))
        inst = DS.CoupledFirstOrderSystem(mesh, 0.5, init_disc=False)
        assert_equal(inst.discs.keys(), tuple())
        inst.initDiscretiser('E', BCs.E)
        self.assert_(isinstance(inst.discs.E, Discretiser.PformDiscretiser))

    def test_reset_history(self):
        self.inst.addLogger('E', [1,2])
        self.inst.addReconstructedLogger('B', [1], N.array([[0,0,0,1]]))
        dofs = self.inst.dofs
        dofs.E.dofArray[:]=1
        dofs.H.dofArray[:]=2
        dofs.D.dofArray[:]=3
        dofs.B.dofArray[:]=4
        self.inst.n=5
        self.inst.reset_history()
        for discname in ('E', 'H', 'D', 'B'):
            self.assert_(N.alltrue(dofs[discname].dofArray == 0))
        assert_equal(self.inst.loggedDOFs.values(),[{}, {}, {}, {}])
        assert_equal(self.inst.loggedReconstructed.values(),[[], [], [], []])


