# Authors:
# Neilen Marais <nmarais@gmail.com>

from __future__ import division

import pickle
import numpy as np
import sys
import dolfin
sys.path.append('../../')
from FenicsCode.Consts import eps0, mu0, c0, Z0
from FenicsCode.Utilities.Converters import as_dolfin_vector
import FenicsCode.Utilities.Optimization

from FenicsCode.Sources.PostProcess import ComplexVoltageAlongLine

# Enable dolfin's form optimizations
FenicsCode.Utilities.Optimization.set_dolfin_optimisation()


fname = 'data/f-1000000000.000000_o-2_s-0.299792_l-0.001000_h-0.166667'
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

complex_voltage = ComplexVoltageAlongLine(V)
complex_voltage.set_dofs(x)
volts = complex_voltage.calculate_voltage(*data['source_endpoints'])

print 'source power: ', volts*data['I'], ' energy flux: ', energy_flux
