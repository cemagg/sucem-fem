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
from NewCode.PostProc import LocatePoints
g_eps = 1e-10                           # Geometrical tollerance
g_contr = 1.-g_eps
x_hat = N.array([1,0,0], N.float64)
y_hat = N.array([0,1,0], N.float64)
z_hat = N.array([0,0,1], N.float64)


def get_ffts(ts_dofs, n, norm=1., fft_range=None):
    no_dofs = ts_dofs.shape[1]
    if fft_range is not None: n_low, n_high = fft_range
    else: n_low = 0 ; n_high = n
    no_fftpts = n_high-n_low
    output = N.zeros((no_fftpts, no_dofs), dtype=N.complex128)
    for i, dof_ts in enumerate(ts_dofs.T):
        output[:,i] = (scipy.fft(dof_ts, n=n)/norm)[n_low:n_high]
    return output

base_filename = '+pml_brick_dipole_FF_h8.5_o2_pml10_dt75.0'

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
E_dofs = E_disc.newDOFs()
B_dofs = B_disc.newDOFs()



surfs = [Struct(const=2, constval=coll_corners[1][2], n_hat=z_hat ),
         Struct(const=2, constval=coll_corners[0][2], n_hat=-z_hat),         
         Struct(const=1, constval=coll_corners[1][1], n_hat=y_hat ),
         Struct(const=1, constval=coll_corners[0][1], n_hat=-y_hat),         
         Struct(const=0, constval=coll_corners[1][0], n_hat=x_hat ),
         Struct(const=0, constval=coll_corners[0][0], n_hat=-x_hat)]
                       
                       

on_collp = on_box_surf(coll_corners[0], coll_corners[1], g_eps)
on_intgsurf_ps = [close_to_point(s.constval, g_eps) for s in surfs]
coll_fs_gen = lambda constcoord, const_p, face: on_collp(face.nodeCoords) and \
          N.all(const_p(face.nodeCoords[:,constcoord]))
coll_fs_s = [partial(coll_fs_gen, surfs[i].const, on_intgsurf_ps[i])
             for i in range(len(surfs))]
in_coll = in_box(*coll_corners)

submeshes = [] ; subGeomEntities_s = []

for coll_fs in coll_fs_s:
    surf_matcher = BrickSurfaceFieldMatcher()
    surf_matcher.initSubdim(E_disc, coll_fs, allfree, dtype=N.complex128)
    subGeomEntities_s.append(surf_matcher.subGeomEntities)
    submeshes.append(surf_matcher.subSurf)

subdiscs = [] ; subdofs = []
for subGeomEntities, submesh in zip(subGeomEntities_s, submeshes):
    subdiscs.append( BrickSubDimDiscretiser.PformOutwardSubDimDiscretiser(
        1, submesh, subGeomEntities, Permuter, E_disc))
    subdiscs[-1].set_interior_selector(in_coll)
    subdiscs[-1].setIntegrationRule(1)
    subdofs.append(subdiscs[-1].newDOFs())

def surf_vals_eval(subdisc, super_dofnos, dof_ffts):
    for dof_fft_n in zip(fft_dof_ns):
        for el in subdisc.elements:
            pv = el.physVals()[:,0]
            perm_l, perm_g = el.permutation()
            eldofs = dof_ffts[dof_fft_n][super_dofnos][perm_g]
            yield N.add.reduce(pv*eldofs[:, N.newaxis])

E_surf_vals = N.array([N.array(
    [val for val in surf_vals_eval(subdisc, super_dofnos, E_dof_ffts)],
    N.complex128).reshape((len(ks), len(subdisc.elements), 3))
             for subdisc, super_dofnos in zip(subdiscs, super_dofnos_s)])

E_surf_vals_max = N.array([[[N.cross(surf.n_hat, pt_val) for pt_val in k_vals]
                            for k_vals in surf_vals]
                           for surf, surf_vals in izip(surfs, dipole_ntff_maxres.surf_E_vals)],
                          N.complex128)
E_surf_vals_max_norm = E_surf_vals_max.copy()
E_surf_vals_max_norm[E_surf_vals_max_norm==0] = 1


E_logdofs =  read_filelog(file(base_filename+'_E.numpyraw', 'rb'))
#B_logdofs =  read_filelog(file(base_filename+'_B.numpyraw', 'rb'))

E_globdoftables_s = [Struct(
    (ent, E_disc.permuter.globalEntityPermutationTable(ent))
    for ent in subdisc.geomEntities)
                     for subdisc in subdiscs]


super_dofnos_s = [N.hstack([E_globdoftables[ent][subdisc.geomEntities[ent][:].superNo].flatten()
                            for ent in sorted(subdisc.geomEntities.keys())])
                  for subdisc in subdiscs]

n_fft = 2**13
dt = run_info.dt ; drv_ts = run_info.drv_ts
f_nyq = 1/2/dt
df = 2*f_nyq/n_fft
freq = N.arange(n_fft)*df
drv_fft = scipy.fft(drv_ts, n=n_fft)
f_min, f_max = 0.25, 1
fft_min, fft_max = int(N.floor(f_min/df)), int(N.ceil(f_max/df))+1
f_mid = .5
fft_mid = int(N.floor(f_mid/df))
fft_ns = N.array([fft_min, fft_mid, fft_max-1], N.int32)
fft_dof_ns = fft_ns - fft_min
freqs = fft_ns*df
ks = freqs*2*N.pi
all_freqs = N.arange(n_fft)*df
E_dof_ffts = get_ffts(E_logdofs, n_fft, norm=drv_fft, fft_range=(fft_min,fft_max))

E_integrals = [E_intg(subdisc) for subdisc in subdiscs]

r_hat = x_hat
Ls = N.array([-N.array([E_integral.calc(r_hat, k, E_fft[super_dofnos])
                for k, E_fft in izip(ks, E_dof_ffts[fft_dof_ns])])
      for E_integral, super_dofnos in zip(E_integrals, super_dofnos_s)])


import dipole_ntff_maxres
from dipole_ntff_maxres import L_surfs

test_pts = N.array([[-1,0,0], [1,0,0],
                    [-1,-0.94117647, -0.94117647],
                    [1,-0.94117647, -0.94117647]], N.float64)*g_contr
elnos, el_coords = LocatePoints(brick_mesh, test_pts)
E_reconparms = list(E_dofs.calc_reconparms(elnos, el_coords))
E_recons = N.array([E_dofs.recon_fromparms(E_reconparms, ed) for ed in E_logdofs])
E_recons_fftdofs = N.array([E_dofs.recon_fromparms(E_reconparms, ed)
                            for ed in E_dof_ffts[fft_dof_ns]])

# #B_reconparms = list(B_dofs.calc_reconparms(elnos, el_coords))
# #B_recons = N.array([B_dofs.recon_fromparms(B_reconparms, ed) for ed in B_logdofs])
# ts = N.arange(len(E_recons))*dt

# from dipole_ntff_maxres import E_z_maxima, H_vals_maxima, E_z_phasor, E_z_phasor_freqs

# E_z_fft = scipy.fft(E_recons[:,0,2], n=n_fft)/drv_fft
# E_z_max_fft = scipy.fft(E_z_maxima, n=n_fft)/drv_fft


# sub_elno = [i for i, sn in enumerate(subdisc.elements[:].superFacenos)
#           if sn in brick_mesh.elements[elnos[0]].facenos][0]
# sub_elcoord = N.array([0.5,0.5], N.float64)

# E_sub_physvals = subdisc.elements[sub_elno].physValsAtPoints([sub_elcoord])[:,0]
# E_sub_perm = subdisc.elements[sub_elno].permutation()[1]
# E_sub_recons = N.array([N.sum(
#     E_sub_physvals*E_fftd[super_dofnos][E_sub_perm][:, N.newaxis], axis=0)
#                         for E_fftd in E_dof_ffts], N.complex128)


# # plot(all_freqs[fft_min:fft_max], N.abs(E_z_fft[fft_min:fft_max]),
# #      label='abs(numeric)')
# plot(all_freqs[fft_min:fft_max], N.abs(E_recons_fftdofs[:,0,2]),
#      label='abs(numeric_fftdofs)')
# plot(all_freqs[fft_min:fft_max], N.abs(E_sub_recons[:,1]),
#      label='abs(numeric_sub)')
# # plot(all_freqs[fft_min:fft_max], N.abs(E_z_max_fft[fft_min:fft_max]),
# #      label='abs(analytical)')
# plot(E_phasor_freqs, N.abs(E_phasor[:,2]), label='abs(anl_phasor)')

# # plot(all_freqs[fft_min:fft_max], N.angle(E_z_fft[fft_min:fft_max]),
# #      label='ang(numeric)')
# plot(all_freqs[fft_min:fft_max], N.angle(E_recons_fftdofs[:,0,2]),
#      label='ang(numeric_fftdofs)')
# plot(all_freqs[fft_min:fft_max], N.angle(E_sub_recons[:,1]),
#      label='ang(numeric_sub)')
# # plot(all_freqs[fft_min:fft_max], N.angle(E_z_max_fft[fft_min:fft_max]),
# #      label='ang(analytical)')
# plot(E_phasor_freqs, N.angle(E_phasor[:,2]), label='ang(anl_phasor)')
# legend(loc='best')
