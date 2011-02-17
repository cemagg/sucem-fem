"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
import pickle
from numpy import *
import numpy as N
import scipy
#
# Local Imports
#
import NewCode
import NewCode.eMAGUSImport as eMAGUSImport
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode import DifferentialForm, Utilities
from NewCode.Utilities import Struct
from NewCode.DifferentialForm import BrickDiscretiser

h0 = 1/5.
a,b,c = 29,23,19

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [a*h0, b*h0, c*h0]))

print 'Mesh elements: ', len(mesh.elements)

class CoupledFirstOrderSystem(NewCode.DiscretisedSystem.CoupledFirstOrderSystem):
    DiscretiserModule = BrickDiscretiser

orders = {'E':(1,True), 'H':(1,True), 'D':(1,True), 'B':(1,True)}
coupledSystem = CoupledFirstOrderSystem(
    mesh, BCs=Struct(E=cb, H=free, D=free, B=cb),
    disc_orders=orders)

E,B = [coupledSystem.discs[dn] for dn in ('E', 'B')]

totalDOFs_E = E.totalDOFs

from NewCode.Integration import BrickTrapzIntegrator

E.set_integrator(BrickTrapzIntegrator)
E.setIntegrationRule(0)
E.D().set_integrator(BrickTrapzIntegrator)
E.D().setIntegrationRule(0)
B.set_integrator(BrickTrapzIntegrator)
B.setIntegrationRule(0)
print "M_e"
M_e = E.matrix.mass()
print "M_b"
M_b = B.matrix.mass()
print "Inverting M_b"
M_b_inv = N.mat(scipy.linalg.inv(M_b.todense(), overwrite_a=True))
print "P_curleb"
P_culreb = E.D().matrix.projectionOnto(B)
P_mu = P_culreb.T
print "C"
C = E.matrix.exteriorDerivative(B)

S = E.matrix.stiffness()
S_P_mu_C = P_culreb.T.matmat(C)
S_P_mu_M_b = P_mu.matmat(M_b_inv*P_mu.T)
