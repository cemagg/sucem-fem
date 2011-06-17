from __future__ import division

import numpy as np
import pickle
import integrate_analytical2 as IA
reload(IA)
IA.cube_size = 0.05
phi = 0
thetas = np.deg2rad(np.linspace(0,90,46))

freq = 1e9
lam = IA.c0/freq
k0 = 2*np.pi*freq/IA.c0
I = 1
l = lam/1000

Ns = []
Ls = []
theta_hats = []
phi_hats = []
E_thetas = []
E_phis = []

for theta in thetas:
    print theta
    r_hat = IA.calc_r_hat(1, theta, phi)      # r_hat doesn't actually depend on r
    theta_hat = IA.calc_theta_hat(1, theta, phi)
    phi_hat = IA.calc_phi_hat(1, theta, phi)
    
    eval_factory = IA.EvalFactory(l, I, k0, r_hat)
    face_integrand = IA.FaceIntegrand(eval_factory)

    N1_integrand = face_integrand.get_face_function('N', 0)
    L1_integrand = face_integrand.get_face_function('L', 0)
    N2_integrand = face_integrand.get_face_function('N', 1)
    L2_integrand = face_integrand.get_face_function('L', 1)
    N3_integrand = face_integrand.get_face_function('N', 2)
    L3_integrand = face_integrand.get_face_function('L', 2)
    N4_integrand = face_integrand.get_face_function('N', 3)
    L4_integrand = face_integrand.get_face_function('L', 3)
    N5_integrand = face_integrand.get_face_function('N', 4)
    L5_integrand = face_integrand.get_face_function('L', 4)
    N6_integrand = face_integrand.get_face_function('N', 5)
    L6_integrand = face_integrand.get_face_function('L', 5)

    ip = IA.IntegrateParts()
    N1 = ip.integrate_all_parts(N1_integrand)
    N2 = ip.integrate_all_parts(N2_integrand)
    N3 = ip.integrate_all_parts(N3_integrand)
    N4 = ip.integrate_all_parts(N4_integrand)
    N5 = ip.integrate_all_parts(N5_integrand)
    N6 = ip.integrate_all_parts(N6_integrand)
    L1 = ip.integrate_all_parts(L1_integrand)
    L2 = ip.integrate_all_parts(L2_integrand)
    L3 = ip.integrate_all_parts(L3_integrand)
    L4 = ip.integrate_all_parts(L4_integrand)
    L5 = ip.integrate_all_parts(L5_integrand)
    L6 = ip.integrate_all_parts(L6_integrand)

    N = N1+N2+N3+N4+N5+N6
    L = L1+L2+L3+L4+L5+L6
    Ns.append(N)
    Ls.append(L)
    theta_hats.append(theta_hat)
    phi_hats.append(phi_hat)

    E_theta_r1 = np.dot(phi_hat, L) + IA.Z0*np.dot(theta_hat, N)
    E_phi_r1 = np.dot(theta_hat, L) - IA.Z0*np.dot(phi_hat, N)
    E_thetas.append(E_theta_r1)
    E_phis.append(E_phi_r1)

pickle.dump(dict(E_thetas=E_thetas, E_phis=E_phis, Ns=Ns, Ls=Ls, theta_hats=theta_hats, phi_hats=phi_hats, thetas=thetas), open('numpy_analytical_integration0.05.pickle', 'w'))

import pylab
pylab.plot(np.rad2deg(thetas), np.abs(E_thetas), label='E_theta')
pylab.plot(np.rad2deg(thetas), np.abs(E_phis), label='E_phis')
pylab.legend(loc='best')
pylab.grid(1)
pylab.show()
