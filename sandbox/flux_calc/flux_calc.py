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
import numpy as np
import sys
import dolfin
sys.path.append('../../')
from sucemfem.Consts import eps0, mu0, c0, Z0
from sucemfem.Utilities.Converters import as_dolfin_vector
import sucemfem.Utilities.Optimization

from sucemfem import Geometry 
from sucemfem.Sources.PostProcess import ComplexVoltageAlongLine
from sucemfem.PostProcessing import CalcEMFunctional

# Enable dolfin's form optimizations
sucemfem.Utilities.Optimization.set_dolfin_optimisation()


fname = 'data/f-1000000000.000000_o-2_s-0.299792_l-0.100000_h-0.166667'
#fname = 'data/f-1000000000.000000_o-2_s-0.299792_l-0.100000_h-0.083333'
data = pickle.load(open(fname+'.pickle'))
mesh = dolfin.Mesh(data['meshfile'])
material_meshfn = dolfin.MeshFunction('uint', mesh, data['materialsfile'])
V = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", data['order'])
x = data['x']
x_r = as_dolfin_vector(x.real)
x_i = as_dolfin_vector(x.imag)
E_r = dolfin.Function(V, x_r)
E_i = dolfin.Function(V, x_i)
k0 = 2*np.pi*data['freq']/c0

n = V.cell().n

ReS = (1/k0/Z0)*dolfin.dot(n, (dolfin.cross(E_r, -dolfin.curl(E_i)) +
                               dolfin.cross(E_i, dolfin.curl(E_r))))*dolfin.ds

energy_flux = dolfin.assemble(ReS)

def boundary(x, on_boundary):
    return on_boundary
E_r_dirich = dolfin.DirichletBC(V, E_r, boundary)
x_r_dirich = as_dolfin_vector(np.zeros(len(x)))
E_r_dirich.apply(x_r_dirich)
E_i_dirich = dolfin.DirichletBC(V, E_i, boundary)
x_i_dirich = as_dolfin_vector(np.zeros(len(x)))
E_i_dirich.apply(x_i_dirich)
x_dirich = x_r_dirich.array() + 1j*x_i_dirich.array()

emfunc = CalcEMFunctional(V)
emfunc.set_E_dofs(x)
emfunc.set_g_dofs(1j*x_dirich.conjugate()/k0/Z0)
emfunc.set_k0(k0)
cell_domains = dolfin.CellFunction('uint', mesh)
cell_domains.set_all(1)
cell_region = 0
boundary_cells = Geometry.BoundaryEdgeCells(mesh)
boundary_cells.mark(cell_domains, cell_region)
emfunc.set_cell_domains(cell_domains, cell_region)

var_energy_flux = emfunc.calc_functional().conjugate()


complex_voltage = ComplexVoltageAlongLine(V)
complex_voltage.set_dofs(x)

volts = complex_voltage.calculate_voltage(*data['source_endpoints'])



print 'source power: ', volts*data['I']
print 'energy flux:      ', energy_flux
print 'var energy flux: ', var_energy_flux
