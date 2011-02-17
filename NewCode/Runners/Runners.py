from __future__ import division
import numpy as N

from NewCode.Utilities import Struct

class TestRun(object):
    drv_fun = None                      # drv_fun(dt, n) -> excitation time waveform
    useLU = True                        # False for iterative solution
    SystemClass = None                  # Class from DiscretisedSystem
    log_divisor = 1
    def __init__(self, mesh, **kwargs):
        self.system = self.SystemClass(mesh, **kwargs)
        self.driveDOFs = self.make_driveDOFs()
        self.system.setDriveDOFs(self.driveDOFs.drive_dofnos,
                                 self.driveDOFs.weights,
                                 self.drv_fun)

    def set_dt(self, dt):
        self.dt = dt
        self.resetHistory()
        if self.useLU: self.system.useLU = True
        print "Setting timestep"
        self.system.setTimestep(dt)
        print "Done"

    def get_stepsCompleted(self):
        return self.system.n

    def resetHistory(self):
        self.system.reset_history()
        self.setupLogging()
        self.system.log()                     # To add an entry for n=0
        
    def make_driveDOFs(self):
        print "Getting source DOFs"
        weights, elPerm = self.system.dofs.calcProjPointfunRHS_with_elPerm(
            matchfun=lambda r: N.array([1,1,1.]), r0=N.array([7.,2.,4.]))
        drive_dofnos = elPerm[1]
        weights = weights[elPerm[0]]
        return Struct(weights=weights, drive_dofnos=drive_dofnos)

    def setupLogging(self):
        divisor = self.log_divisor
        # setup reconstructed or DOF loggers here, eg. calls to
        # self.system.addReconstructedLogger or
        # self.system.addLogger
        pass
    
    def runSteps(self, n):
        self.system.step(n)

    def getResult(self):
        # Extract the results from the logged stuff in the format you want it
        # to be saved by getResults
        pass    

    def yieldResults(self, dt_times):
            for dt, times in dt_times.iteritems():
                self.set_dt(dt)
                for time in sorted(times):
                    prev_steps = self.get_stepsCompleted()
                    self.runSteps(int(time/dt)-prev_steps)
                    yield (dt, time), self.getResult()

    def getResults(self, dt_times):
        return dict((k,v) for k,v in self.yieldResults(dt_times))
