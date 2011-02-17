import numpy as N

from NewCode import PostProc, DifferentialForm, Waveforms

from wg_dof_matching_stuff import RoughCoaxMatcher

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
dt = 0.01

drv_fun = Waveforms.get_d_gaussian(fc=0.5, tpr=-50)

def setupMatcher(mesh):
    global wgm, fs, onPort, direction, xypos, test_pts, test_elnos, test_el_coords
    wgm = RoughCoaxMatcher(mesh)
    fs = wgm.fs
    onPort = wgm.onPort
    direction = N.array([0, 0, 1.], N.float64)
    xypos = [(wgm.r_inner + wgm.r_outter)/2, 0, 0]
    test_pts = PostProc.MakeLine(0.25*direction, 9.75*direction, 401) + xypos
    test_elnos, test_el_coords = PostProc.LocatePoints(mesh, test_pts)
    
