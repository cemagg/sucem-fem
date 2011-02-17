from __future__ import division

from itertools import izip

import pickle
import scipy
import numpy as N

from NewCode.Utilities import partial, close_to_point, Struct, \
     in_box, on_box_surf
from NewCode.IOUtils import read_filelog
from NewCode.DifferentialForm import BrickDiscretiser
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.BrickSubDimMesh import SubSurface
from NewCode.DifferentialForm import BrickSubDimDiscretiser, \
     BrickSubDimDiscretiserEntities, allfree
from NewCode.DifferentialForm.Discretiser import Permuter
from NewCode.Feeds import BrickSurfaceFieldMatcher
from NewCode.NTFF import E_intg
g_eps = 1e-10                           # Geometrical tollerance

def get_ffts(ts_dofs, n, norm=1., fft_range=None):
    no_dofs = ts_dofs.shape[1]
    if fft_range is not None: n_low, n_high = fft_range
    else: n_low = 0 ; n_high = n
    no_fftpts = n_high-n_low
    output = N.zeros((no_fftpts, no_dofs), dtype=N.complex128)
    for i, dof_ts in enumerate(ts_dofs.T):
        output[:,i] = (scipy.fft(dof_ts, n=n)/norm)[n_low:n_high]
    return output

base_filename = '+pml_brick_cube_RCS_h8.0_o1_pml5_dt75.0'

run_info = pickle.load(file(base_filename+'.pickle'))

brick_mesh = BrickMesh.Mesh(run_info.brick_listmesh)
print 'Brick_mesh elements: ', len(brick_mesh.elements)
order = run_info.order
coll_corners = run_info.Gamma_C_corners

free_edges = set(run_info.ntff.log_ents.edge)
free_faces = set(run_info.ntff.log_ents.face)
free_vols = set(run_info.ntff.log_ents.vol)

def log_free(ent):
    if len(ent.nodes) == 2: return ent.index in free_edges
    if len(ent.nodes) == 4: return ent.index in free_faces
    raise Exception('oops, seems not to be a brick edge/face')

def log_volfree(vol):
    return vol.index in free_vols

dead_elements = set(xrange(len(brick_mesh.elements))) - free_vols

E_disc = BrickDiscretiser.setup_PformDiscretiser(
    brick_mesh, 1, order=order, freeFun=log_free,
    vol_freeFun=log_volfree, btype='cohen98', dead_elements=dead_elements)
B_disc = BrickDiscretiser.setup_PformDiscretiser(
    brick_mesh, 2, order=order, freeFun=log_free,
    vol_freeFun=log_volfree, btype='cohen98', dead_elements=dead_elements)
E_disc.diagonalise()
B_disc.diagonalise()

on_collp = on_box_surf(coll_corners[0], coll_corners[1], g_eps)
coll_fs = lambda face: on_collp(face.nodeCoords)
in_coll = in_box(*coll_corners)
# submesh = SubSurface(brick_mesh, coll_fs)

# sdEntityList = BrickSubDimDiscretiserEntities.SubDimDiscretiserEntityList
# sdEdge = BrickSubDimDiscretiserEntities.Edge
# sdFace = BrickSubDimDiscretiserEntities.Face

surf_matcher = BrickSurfaceFieldMatcher()
surf_matcher.initSubdim(E_disc, coll_fs, allfree, dtype=N.complex128)
submesh = surf_matcher.subSurf
subGeomEntities = surf_matcher.subGeomEntities
# subGeomEntities = Struct(edge=sdEntityList(sdEdge(mesh=submesh, freefun=allfree,
#                                                   attrs=submesh.edges.list_repr())))

# if 'face' in E_disc.basisSet.fns:
#     subGeomEntities['face'] = sdEntityList(sdFace(
#         mesh=submesh, freefun=allfree, attrs=submesh.elements.list_repr()))

subdisc = BrickSubDimDiscretiser.PformOutwardSubDimDiscretiser(
    1, submesh, subGeomEntities, Permuter, E_disc)
subdisc.set_interior_selector(in_coll)
subdofs = subdisc.newDOFs()




E_logdofs =  read_filelog(file(base_filename+'_E.numpyraw', 'rb'))
B_logdofs =  read_filelog(file(base_filename+'_B.numpyraw', 'rb'))

E_globdoftables = Struct(
    (ent, E_disc.permuter.globalEntityPermutationTable(ent))
    for ent in subdisc.geomEntities)


super_dofnos = N.hstack([E_globdoftables[ent][subdisc.geomEntities[ent][:].superNo].flatten()
                         for ent in sorted(subdisc.geomEntities.keys())])

freq_step = 1
n_fft = 2**13
dt = run_info.dt ; drv_ts = run_info.drv_ts
f_nyq = 1/2/dt
df = 2*f_nyq/n_fft
freq = N.arange(n_fft)*df
drv_fft = scipy.fft(drv_ts, n=n_fft)
f_min, f_max = 0.4, 2
fft_min, fft_max = int(N.floor(f_min/df)), int(N.ceil(f_max/df))+1

E_dof_ffts = get_ffts(E_logdofs, n_fft, norm=drv_fft, fft_range=(fft_min,fft_max))
B_dof_ffts = get_ffts(B_logdofs, n_fft, norm=drv_fft, fft_range=(fft_min,fft_max))

E_integral = E_intg(subdisc)
ks = 2*N.pi*N.arange(fft_min, fft_max)*df
freqs = N.arange(fft_min, fft_max)[0::freq_step]*df
r_hat = N.array([0,0,-1.])
vals = N.array([N.cross(r_hat, E_integral.calc(r_hat, k, E_fft[super_dofnos]))
                for k, E_fft in izip(ks[0::freq_step], E_dof_ffts[0::freq_step])])
M_e = E_disc.matrix.mass()
M_b = B_disc.matrix.mass()
C = E_disc.matrix.partialExteriorDerivative(B_disc)

x_hat = N.array([1.,0,0])
y_hat = N.array([0,1.,0])
u_hats = [x_hat, y_hat]

def matchfn(u_hat, k, r_hat, r):
    return u_hat*N.exp(1j*k*N.dot(r, r_hat))

def a_E(e,b,v,w):
    return 1j*w*N.dot(v, M_e.matvec(e)) + N.dot(v, C.T.matvec(M_b.matvec(b)))

def a_E_dt1(e,b,v,w):
    return (N.exp(-1j*dt/2) - N.exp(1j*dt/2))/dt*N.dot(v, M_e.matvec(e)) \
           + N.dot(v, C.T.matvec(M_b.matvec(b)))

def a_E_dt2(e,b,v,w):
    return (N.exp(-1j*dt/2) - N.exp(1j*dt/2))/dt*N.dot(v, M_e.matvec(e)) \
           - N.dot(v, C.T.matvec(M_b.matvec(b)))

H_vals = N.zeros_like(vals)
H_vals_dt1 = N.zeros_like(vals)
H_vals_dt2 = N.zeros_like(vals)
subM_solve = scipy.sparse.linalg.dsolve.factorized(
    surf_matcher.DOFs.matrix.mass().astype(N.complex128))
H_RHSes = N.array([surf_matcher.DOFs.calcProjRHS(
    partial(matchfn, u_hat, k, r_hat))
                   for k in ks[0::freq_step] for u_hat in u_hats])


for i, (k, E_fft, B_fft) in enumerate(
    izip(ks[0::freq_step], E_dof_ffts[0::freq_step], B_dof_ffts[0::freq_step])):
    print i
    for u_j, u_hat in enumerate(u_hats):
        v = N.zeros_like(E_fft)
        #         matchfn_k = partial(matchfn, u_hat, k, r_hat)
        #         RHS = surf_matcher.DOFs.calcProjRHS(matchfn_k)
        RHS = H_RHSes[i*len(u_hats)+u_j]
        v[super_dofnos] = subM_solve(RHS)
        H_vals[i]     += u_hat*a_E    (E_fft, B_fft, v, k)
        H_vals_dt1[i] += u_hat*a_E_dt1(E_fft, B_fft, v, k)
        H_vals_dt2[i] += u_hat*a_E_dt2(E_fft, B_fft, v, k)

from NewCode.Consts import c0

H_vals_sp = N.exp(1j*dt/2*ks[0::freq_step])[:,N.newaxis]*H_vals
H_vals_sn = N.exp(-1j*dt/2*ks[0::freq_step])[:,N.newaxis]*H_vals
H_vals_dt1sp = N.exp(1j*dt/2*ks[0::freq_step])[:,N.newaxis]*H_vals_dt1
H_vals_dt1sn = N.exp(-1j*dt/2*ks[0::freq_step])[:,N.newaxis]*H_vals_dt1
H_vals_dt2sp = N.exp(1j*dt/2*ks[0::freq_step])[:,N.newaxis]*H_vals_dt2
H_vals_dt2sn = N.exp(-1j*dt/2*ks[0::freq_step])[:,N.newaxis]*H_vals_dt2

from pylab import plot
#plot(freqs*c0, (ks**2/N.pi/4)*N.abs(vals[:,0])**2)
