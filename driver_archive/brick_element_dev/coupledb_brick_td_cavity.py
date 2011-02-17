"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
import pickle
import numpy as N
import scipy
import random
from scipy import signal
#
# Local Imports
#
import NewCode
from NewCode.Utilities import Struct, partial
from NewCode.Meshes import BrickMesh
from NewCode import DifferentialForm, Waveforms
from NewCode.DifferentialForm import BrickDiscretiser, DiscretiserDOFs
from NewCode.DiscretisedSystem import BrickCoupledFirstOrderSystemB
import brick_cavity

mesh = BrickMesh.Mesh(
    brick_cavity.make_rect_cavity_brick_listmesh(1,0.75,0.5,1/9.))

print 'Mesh elements: ', len(mesh.elements)

weights = 1.0
dt = 0.01
cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
#orders = {'E':(2,True), 'H':(2,True), 'D':(1,False), 'B':(1,False)}
# All lowest order
orders = {'E':(1,True), 'B':(1,True)}
coupledSystem = BrickCoupledFirstOrderSystemB(
    mesh, dt=dt, BCs=Struct(E=cb, B=cb), disc_orders=orders)

totalDOFs_E = coupledSystem.discs.E.totalDOFs
workDOFs = min(int(N.ceil(totalDOFs_E/100.)), 1000)
# drive_dofnos and test_dofnos are randomly chosen
print "Choosing drive/log DOFs"
dofnolist = range(coupledSystem.discs.E.totalDOFs)
random.shuffle(dofnolist)
dofnolist = dofnolist[0:2*workDOFs]
drive_dofnos = dofnolist[0:workDOFs]
test_dofnos = dofnolist[workDOFs:2*workDOFs]
print "Done choosing"

drv_fun = Waveforms.get_gausspulse(1.25, bw=0.9)
coupledSystem.setDriveDOFs_J(drive_dofnos, weights, drv_fun)
coupledSystem.addLogger('E', test_dofnos)
dofs = coupledSystem.dofs
from NewCode.MatrixUtils import add_diagonal_preconditioner
#add_diagonal_preconditioner(dofs.E.matrix.mass())
DiscretiserDOFs.PformDiscretiserBase.useLU = True
coupledSystem.step(2048*2)

from AnalyticResults import PEC_cavity
an_freqs = N.sqrt(PEC_cavity['rect1x0.75x0.5'])/N.pi/2

E_dofs = N.array(coupledSystem.loggedDOFs.E.values()[0]['vals']).transpose()

import analyse_fftstuff
analyse_this = analyse_fftstuff.AnalyseThingy(E_dofs, dt, drv_fun)

analyse_this.plotFreqbars(an_freqs)
analyse_this.xmin = 120*2
analyse_this.xmax = 265*2
analyse_this.plotRange(20*log10(N.sum(N.abs(analyse_this.dfts_series), axis=0)/N.abs(analyse_this.fs)))

#analyse_this.plotRange(20*log10(N.sqrt(N.sum(N.abs(analyse_this.dfts_series)**2, axis=0))))

dump_res = dict(y=20*log10(N.sum(N.abs(analyse_this.dfts_series), axis=0)/N.abs(analyse_this.fs)
                    )[analyse_this.xmin:analyse_this.xmax],
                x=N.arange(analyse_this.xmin, analyse_this.xmax)*analyse_this.df)

#pickle.dump(dump_res, file('+coupledb_brick_td_cavity-' + mesh.FemmeshFilename + '.pickle', 'w'))
