## Copyright (C) 2011 Stellenbosch University
##
## This file is part of SUCEM.
##
## SUCEM is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SUCEM is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with SUCEM. If not, see <http://www.gnu.org/licenses/>. 
##
## Contact: cemagga@gmail.com 
# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import pickle
import numpy as N
#import dolfin
from dolfin import *

import sys
sys.path.append('../../')
from sucemfem.Consts import Z0, c0
from sucemfem.Utilities.MeshGenerators import get_centred_cube
import sucemfem.Utilities.Optimization
from sucemfem.PostProcessing import surface_ntff, variational_ntff

# Enable dolfin's form optimizations
sucemfem.Utilities.Optimization.set_dolfin_optimisation()

# Near-field of an infintesimal dipole
fname = 'reference_dofs-2-0.149896229-0.0499654096667.pickle'
theta_deg = N.linspace(0, 180, 181)
no_ff_pts = len(theta_deg)
phi_deg = N.zeros(no_ff_pts)
data = pickle.load(open(fname))
print data

lam = c0/data['freq']
k0 = data['freq']*2*N.pi/c0
mesh = get_centred_cube(data['domain_size'], data['max_edge_len'])
#plot(mesh,interactive=True)
V = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", data['order'])

#------------------------------
# Calculate using the standard surface integral NTFF
surf_ntff = surface_ntff.NTFF(V)
surf_ntff.set_dofs(data['x'])
surf_ntff.set_frequency(data['freq'])
surf_E_ff = N.array([surf_ntff.calc_pt(th_deg, ph_deg)
                for th_deg, ph_deg in zip(theta_deg, phi_deg)])
surf_E_theta = surf_E_ff[:,0]
surf_E_phi = surf_E_ff[:,1]
#------------------------------
# Calculate using the variational integral NTFF
var_ntff = variational_ntff.NTFF(V)
var_ntff.set_k0(k0)
var_ntff.set_dofs(data['x'])
var_ntff.set_frequency(data['freq'])
var_E_ff = N.array([var_ntff.calc_pt(th_deg, ph_deg)
                for th_deg, ph_deg in zip(theta_deg, phi_deg)])
var_E_theta = var_E_ff[:,0]
var_E_phi = var_E_ff[:,1]
#------------------------------
# Analytical solution of the far-field of an infintesimal dipole
E_theta_an = Z0*N.exp(-1j*k0)*(1j*k0*1*lam/1000)/(4*N.pi)*N.sin(N.deg2rad(theta_deg))
#------------------------------
# Plot results
import matplotlib as mpl
import matplotlib.pyplot as plt
font = {'size'   : 18}

mpl.rc('font', **font)

plt.figure()
plt.plot(theta_deg, N.abs(surf_E_theta), 'ko',markersize=10,markerfacecolor='None')
#plt.plot(theta_deg, N.abs(surf_E_phi), label='surf_ntff |E_phi|')
plt.plot(theta_deg, N.abs(var_E_theta), 'kx', linewidth=1.5, markersize=10)
#plt.plot(theta_deg, N.abs(var_E_phi), label='var_ntff |E_phi|')
plt.plot(theta_deg, N.abs(E_theta_an), '-k',linewidth=1.5)
plt.xlabel('Degrees [$^o$]')
plt.title('Surface Integral versus Variational NTFF')
plt.ylabel('$E_{\Theta}$ [V/m]')
plt.legend(['Surface NTFF','Variational NTFF','Analytical'],loc='lower center')
plt.grid(True)
plt.show()

#import pylab 
#pylab.figure()
#pylab.plot(theta_deg, N.abs(surf_E_theta), '-ko', mfc="None", markersize=8, label='Surface NTFF')
##pylab.plot(theta_deg, N.abs(surf_E_phi), label='surf_ntff |E_phi|')
##pylab.plot(theta_deg, N.abs(var_E_theta), '-rx', markersize=8 , label='Variational NTFF')
##pylab.plot(theta_deg, N.abs(var_E_phi), label='var_ntff |E_phi|')
#pylab.plot(theta_deg, N.abs(E_theta_an), '--b', label='Analytical')
#pylab.xlabel('Degrees [$^o$]')
#pylab.title('Surface Integral vs. Variational NTFF')
#pylab.ylabel('$E_{\Theta}$ [V/m]')
#pylab.legend(loc='lower center')
#pylab.grid(1)
#pylab.show()
