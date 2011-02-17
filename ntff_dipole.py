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
freqs = fft_ns*df
E_dof_ffts = get_ffts(E_logdofs, n_fft, norm=drv_fft, fft_range=(fft_min,fft_max)
                      )[fft_ns - fft_min]
B_dof_ffts = get_ffts(B_logdofs, n_fft, norm=drv_fft, fft_range=(fft_min,fft_max)
                      )[fft_ns - fft_min]

E_integral = E_intg(subdisc)
ks = 2*N.pi*freqs

r_hat = N.array([1,0,0.])
print "Calculating E Integrals"
Ls = -N.array([E_integral.calc(r_hat, k, E_fft[super_dofnos])
               for k, E_fft in izip(ks, E_dof_ffts)])
M_e = E_disc.matrix.mass()
M_b = B_disc.matrix.mass()
C = E_disc.matrix.partialExteriorDerivative(B_disc)

z_hat = N.array([0.,0,1])
y_hat = N.array([0,1.,0])
u_hats = [z_hat, y_hat]

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

Ns = N.zeros_like(Ls)
Ns_dt1 = N.zeros_like(Ls)
Ns_dt2 = N.zeros_like(Ls)
subM_solve = scipy.sparse.linalg.dsolve.factorized(
    surf_matcher.DOFs.matrix.mass().astype(N.complex128))
H_RHSes = N.array([surf_matcher.DOFs.calcProjRHS(
    partial(matchfn, u_hat, k, r_hat))
                   for k in ks for u_hat in u_hats])

print "Calculating B integrals"
for i, (k, E_fft, B_fft) in enumerate(
    izip(ks, E_dof_ffts, B_dof_ffts)):
    print i
    for u_j, u_hat in enumerate(u_hats):
        v = N.zeros_like(E_fft)
        #         matchfn_k = partial(matchfn, u_hat, k, r_hat)
        #         RHS = surf_matcher.DOFs.calcProjRHS(matchfn_k)
        RHS = H_RHSes[i*len(u_hats)+u_j]
        v[super_dofnos] = subM_solve(RHS)
        Ns[i]     += u_hat*a_E    (E_fft, B_fft, v, k)
        Ns_dt1[i] += u_hat*a_E_dt1(E_fft, B_fft, v, k)
        Ns_dt2[i] += u_hat*a_E_dt2(E_fft, B_fft, v, k)

from NewCode.Consts import c0

Ns_sp = N.exp(1j*dt/2*ks)[:,N.newaxis]*Ns
Ns_sn = N.exp(-1j*dt/2*ks)[:,N.newaxis]*Ns
Ns_dt1sp = N.exp(1j*dt/2*ks)[:,N.newaxis]*Ns_dt1
Ns_dt1sn = N.exp(-1j*dt/2*ks)[:,N.newaxis]*Ns_dt1
Ns_dt2sp = N.exp(1j*dt/2*ks)[:,N.newaxis]*Ns_dt2
Ns_dt2sn = N.exp(-1j*dt/2*ks)[:,N.newaxis]*Ns_dt2

from pylab import plot
#plot(freqs*c0, (ks**2/N.pi/4)*N.abs(vals[:,0])**2)

### Results from Maxima

Ls_max = N.array([[[0.0,-0.35270705335917*1j-0.65257726581273,-2.720046410331633*10**-19*1j-4.9096529311201388*10**-19],[0.27453599789206*1j+0.45344667336181,-0.27453599789206*1j-0.45344667336181,2.801578413702544*10**-19-4.3526139482786949*10**-19*1j],[0.35270705335917*1j+0.65257726581273,-1.3877787807814457*10**-17*1j,1.2829243840112919*10**-19*1j+1.171902081548777*10**-19],[0.0,-0.28838444740779*1j-0.40119901677968,-1.0177044392397268*10**-18*1j-9.97465998686664*10**-19],[0.28838444740779*1j+0.40119901677968,0.0,1.7347234759768071*10**-18],[0.0,0.0,0.0]],[[0.0,0.32104842827596*1j-0.0010957493662258,0.0],[0.39753907892604-0.39875551136918*1j,0.39875551136918*1j-0.39753907892604,-3.3989015305972237*10**-18*1j-3.2954349126078512*10**-18],[0.0010957493662258-0.32104842827596*1j,-6.9388939039072284*10**-18,0.0],[0.0,0.078076807912481*1j-0.37919521244395,8.3729319773813881*10**-18*1j+3.3306690738754691*10**-18],[0.37919521244395-0.078076807912482*1j,0.0,-3.4694469519536142*10**-18*1j],[6.9388939039072284*10**-18,-6.9388939039072284*10**-18,0.0]],[[-1.3877787807814457*10**-17*1j,0.0082378208502777*1j+0.28672000959199,0.0],[0.26830533829659*1j+0.63952675862079,-0.26830533829659*1j-0.63952675862079,-2.1300091320360551*10**-18*1j-5.575728921779363*10**-18],[-0.0082378208502777*1j-0.28672000959199,6.9388939039072284*10**-18*1j,0.0],[0.0,-0.018212988578969*1j-0.67335590471424,5.551115123125783*10**-18*1j+5.0885221961986346*10**-18],[0.018212988578969*1j+0.67335590471424,-1.3877787807814457*10**-17*1j-1.3877787807814457*10**-17,-3.4694469519536142*10**-18*1j],[-3.4694469519536142*10**-17*1j-1.3877787807814457*10**-17,3.8163916471489756*10**-17*1j+2.7755575615628914*10**-17,0.0]]])
Ls_max[N.abs(Ls_max) < 1e-5] = 0

Ns_max = N.array([[[0.0,0.0,0.34742273879211-0.35270706199386*1j],[0.0,0.0,0.35872960080766-0.3882525019617*1j],[0.0,0.0,0.34742273879211-0.35270706199386*1j],[0.09616965358641-0.14501842240589*1j,-2.9305261920834869*10**-19*1j-7.5348870661369921*10**-18,0.52878850085621-0.55285564693721*1j],[-2.0238440553062755*10**-19*1j-9.8069700508555494*10**-18,0.09616965358641-0.14501842240589*1j,0.52878850085621-0.55285564693721*1j],[-8.1416355139178147*10**-18*1j-6.661338147750939*10**-18,-3.2763915037826836*10**-18*1j-1.1023281053389332*10**-17,0.69484547758422-0.70541412398771*1j]],[[0.0,0.0,0.32104881720144*1j+0.99890456425829],[0.0,0.0,0.56392435682461*1j+0.43779437748422],[0.0,0.0,0.32104881720144*1j+0.99890456425829],[-0.11473591901433*1j-0.34562023562477,5.0599956323714607*10**-18*1j+5.0416461129366822*10**-18,0.11811671614325-0.0043177883071575*1j],[8.5117098554595328*10**-18*1j-3.700743415417193*10**-19,-0.11473591901433*1j-0.34562023562477,0.11811671614325-0.0043177883071576*1j],[0.0,0.0,0.002193405924227-0.64209445136311*1j]],[[0.0,0.0,0.0082355269111796*1j+1.286726577322081],[0.0,0.0,0.095586569034019-0.37944279580089*1j],[0.0,0.0,0.0082355269111796*1j+1.286726577322081],[-0.27518449906288*1j-0.25629986973208,2.9898922843724641*10**-19*1j+6.0435067444291517*10**-18,-0.30094512511054*1j-0.20857208583355],[9.6219328800846903*10**-18*1j-4.440892098500627*10**-18,-0.27518449906288*1j-0.25629986973208,-0.30094512511054*1j-0.20857208583355],[0.0,0.0,-0.016489706529471*1j-0.57339977312137]]])
Ns_max[N.abs(Ns_max) < 1e-5] = 0

E_inf_ans = N.array([[[0,0,-1j/8],[0,0,-1j/8],[0,0,-1j/8],[1j/16,0,-1j/16],[0,1j/16,-1j/16],[0,0,0]],[[0,0,-1j/4],[0,0,-1j/4],[0,0,-1j/4],[1j/8,0,-1j/8],[0,1j/8,-1j/8],[0,0,0]],[[0,0,-1j/2],[0,0,-1j/2],[0,0,-1j/2],[1j/4,0,-1j/4],[0,1j/4,-1j/4],[0,0,0]]])

