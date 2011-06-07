from __future__ import division

import numpy as N
import dolfin
from dolfin import curl, cross, dx, ds, Constant, dot
import pickle

import sys
sys.path.append('../../')
from FenicsCode.Utilities.MeshGenerators import get_centred_cube
from FenicsCode.Utilities.MeshIO import femmesh_2_dolfin_mesh
from FenicsCode.Consts import Z0, c0

theta_deg = N.linspace(0, 180, 91)
phi_deg = 0

#fname = 'dofs-2-0.599584916-0.0499654096667.pickle'
#fname = 'dofs-2-1.199169832-0.0499654096667.pickle'
fname = 'interpdofs-3-0.599584916-0.0499654096667.pickle'
#fname = 'dofs_sph-2-sphere-r1m-6.pickle'
#fname = 'dofs_sph-2lam-2-sphere-r1m-15.pickle'

data = pickle.load(open(fname))
lam = c0/data['freq']
k0 = data['freq']*2*N.pi/c0
# Get real and imaginary parts of the DOF vector
x_r = N.real(data['x']).copy()
x_i = N.imag(data['x']).copy()
# Set up mesh and function space
mesh = get_centred_cube(data['domain_size'], data['max_edge_len'])
# mesh_file = '../solvers/meshes/%s.femmesh' % data['mesh_id']
# mesh = femmesh_2_dolfin_mesh(mesh_file)
# mesh.coordinates()[:] *= 2*lam 
V = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", data['order'])
# Set up real and imaginary E fields
E_r = dolfin.Function(V)
E_i = dolfin.Function(V)
E_r.vector()[:] = x_r
E_i.vector()[:] = x_i

# Outward pointing normal
n = V.cell().n
# \vec{r'}, i.e. rprime is simply the position vector
rprime = V.cell().x

#------------------------------
# Set up eqivalent magnetic current forms
M_r = -cross(n, E_r)
M_i = -cross(n, E_i)
#------------------------------
# Set up magnetic field and equivalent electric current forms
H_r = -curl(E_i)/(k0*Z0)
H_i = curl(E_r)/(k0*Z0)
J_r = cross(n, H_r)
J_i = cross(n, H_i)
#------------------------------
# Calculate various unit vectors
E_theta = []
E_phi = []
L_theta_vals = []
L_phi_vals = []
N_theta_vals = []
N_phi_vals = []
for theta in N.deg2rad(theta_deg):
    phi = N.deg2rad(phi_deg)
    rhat_ = N.array([N.sin(theta)*N.cos(phi),
                     N.sin(theta)*N.sin(phi),
                     N.cos(theta)], dtype=N.float64)
    theta_hat_ = N.array([N.cos(theta)*N.cos(phi),
                          N.cos(theta)*N.sin(phi),
                          -N.sin(theta)], dtype=N.float64)
    phi_hat_ = N.array([-N.sin(phi), N.cos(phi), 0.], dtype=N.float64)
    rhat = Constant(list(rhat_))
    theta_hat = Constant(list(theta_hat_))
    phi_hat = Constant(list(phi_hat_))
    # phase term to be used in sin/cos
    phase = k0*dot(rprime, rhat)
    #------------------------------
    # Set up form for far field potential N
    N_r = J_r*dolfin.cos(phase) - J_i*dolfin.sin(phase)
    N_i = J_r*dolfin.sin(phase) + J_i*dolfin.cos(phase)
    # UFL does not seem to like vector valued functionals, so we split the
    # final functional from into theta and phi components
    N_r_theta = dot(theta_hat, N_r)*ds
    N_r_phi = dot(phi_hat, N_r)*ds
    N_i_theta = dot(theta_hat, N_i)*ds
    N_i_phi = dot(phi_hat, N_i)*ds

    #------------------------------
    # Set up form for far field potential L
    L_r = M_r*dolfin.cos(phase) - M_i*dolfin.sin(phase)
    L_i = M_r*dolfin.sin(phase) + M_i*dolfin.cos(phase)
    L_r_theta = dot(theta_hat, L_r)*ds
    L_r_phi = dot(phi_hat, L_r)*ds
    L_i_theta = dot(theta_hat, L_i)*ds
    L_i_phi = dot(phi_hat, L_i)*ds

    #------------------------------
    # Now evaluate numerically, adding the real and imaginary parts
    L_theta = dolfin.assemble(L_r_theta) + 1j*dolfin.assemble(L_i_theta)
    L_phi = dolfin.assemble(L_r_phi) + 1j*dolfin.assemble(L_i_phi)
    N_theta = dolfin.assemble(N_r_theta) + 1j*dolfin.assemble(N_i_theta)
    N_phi = dolfin.assemble(N_r_phi) + 1j*dolfin.assemble(N_i_phi)    
    N_theta_vals.append(N_theta)
    N_phi_vals.append(N_phi)
    L_theta_vals.append(L_theta)
    L_phi_vals.append(L_phi)

    #------------------------------
    # Calculate the far fields normalised to radius 1.
    r_fac = 1j*k0*N.exp(-1j*k0)/(4*N.pi)
    E_theta.append(-r_fac*(L_phi + Z0*N_theta))
    E_phi.append(r_fac*(L_theta - Z0*N_phi))

E_theta = N.array(E_theta)
E_phi = N.array(E_phi)

pickle.dump(dict(E_theta=E_theta, E_phi=E_phi),
            open(fname.replace('dofs', 'results', 1), 'w'))

import pylab 
pylab.plot(theta_deg, N.abs(E_theta), label='|E_theta|')
pylab.plot(theta_deg, N.abs(E_phi), label='|E_phi|')
E_theta_an = Z0*N.exp(-1j*k0)*(1j*k0*1*lam/1000)/(4*N.pi)*N.sin(N.deg2rad(theta_deg))
pylab.plot(theta_deg, N.abs(E_theta_an), label='analytical')
pylab.legend()
pylab.grid(1)
pylab.show()
