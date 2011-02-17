"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
import pickle
import random
import numpy as N
import scipy
#
# Local Imports
#
import NewCode
from NewCode.Meshes import BrickMesh
from NewCode.Utilities import Struct, partial
from NewCode import DifferentialForm, Waveforms
from NewCode.DifferentialForm import BrickDiscretiser
from NewCode.DiscretisedSystem import BrickCurlCurlNewmark

import brick_cavity
mesh = BrickMesh.Mesh(
    brick_cavity.make_rect_cavity_brick_listmesh(1,0.75,0.5,1/3.))




print 'Mesh elements: ', len(mesh.elements)

weights = 1.0
dt = 0.01
cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
newmarkSystem = BrickCurlCurlNewmark(mesh, order=5, BC=cb)
totalDOFs = newmarkSystem.disc.totalDOFs
workDOFs = min(int(N.ceil(totalDOFs/100.)), 1000)
# drive_dofnos and test_dofnos are randomly chosen
print "Choosing drive/log DOFs"
dofnolist = range(newmarkSystem.disc.totalDOFs)
random.shuffle(dofnolist)
#dofnolist = pickle.load(file('+good_dofnolist.pickle'))
dofnolist = dofnolist[0:2*workDOFs]
drive_dofnos = dofnolist[0:workDOFs]
test_dofnos = dofnolist[workDOFs:2*workDOFs]
print "Done choosing"

drv_fun = Waveforms.get_gausspulse(1.25, bw=0.9)
newmarkSystem.setDriveDOFs(drive_dofnos, weights, drv_fun)
newmarkSystem.addLogger(test_dofnos)
newmarkSystem.useLU = True
newmarkSystem.setTimestep(dt)
newmarkSystem.step(4096)

from AnalyticResults import PEC_cavity
an_freqs = N.sqrt(PEC_cavity['rect1x0.75x0.5'])/N.pi/2

E_dofs = N.array(newmarkSystem.loggedDOFs.values()[0]['vals']).transpose()

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

#pickle.dump(dump_res, file('+newmark_cavity-' + mesh.FemmeshFilename + '.pickle', 'w'))


"""
Plotty stuff:

import analyse_fftstuff
E_dofs = N.array(newmarkSystem.loggedDOFs.values()[0]['vals']).transpose()
analyse_this = analyse_fftstuff.AnalyseThingy(E_dofs, dt, drv_fun)
analyse_this.plotFreqbars(an_freqs)
analyse_this.plotRange(20*log10(N.sum(N.abs(analyse_this.dfts_series), axis=0)/N.abs(analyse_this.fs)))

"""


"""
Results:

Using CT/LN newmark beta=0.3, dt = 0.02
                freq

mesh edge len   0.83333333,1.11803399,1.20185043x2,1.30170828x2,1.41421356,1.42400062,1.56347192x2

0.5/6             1/546       3/733      4/788          4/853      5/927       11/933     6/1024
0.5/12            

"""
