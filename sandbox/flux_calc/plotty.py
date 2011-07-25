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
## along with SUCEM. If not, see <http:##www.gnu.org/licenses/>. 
##
## Contact: cemagga@gmail.com 
from __future__ import division

import pylab
import numpy as np
import run_data
reload(run_data)

def get_unit_gradients(resdict, ind=0, op=lambda x: x):
    return dict((k, np.gradient(op(v[ind])))
                for k,v in resdict.iteritems())

def get_gradients(resdict, op1=lambda x: x, op2=lambda x: x):
    return dict((k, np.gradient(op1(v[1])) / np.gradient(op2(v[0])))
                for k,v in resdict.iteritems())

def logify(resdict, op=lambda x: x):
    return dict((k, np.log(op(v))) for k,v in resdict.iteritems())

vflux_gradients_unity = get_unit_gradients(run_data.vflux, ind=1,
                                           op=lambda x: np.real(x))
vflux_gradients = get_gradients(run_data.vflux, op1=np.real)
vflux_log = logify(vflux_gradients, op=np.abs)
vflux_log_h = logify(run_data.vflux, op=lambda x: x[0])

pylab.figure(1)
pylab.hold(0)
pylab.plot(-vflux_log_h['1r'], vflux_log['1r'], label='reversed')
pylab.hold(1)
pylab.plot(-vflux_log_h['1'], vflux_log['1'], label='forward')
pylab.legend(loc=0)
pylab.figure(2)
pylab.hold(0)
pylab.plot(-vflux_log_h['2r'], vflux_log['2r'], label='reversed')
pylab.hold(1)
pylab.plot(-vflux_log_h['2'], vflux_log['2'], label='forward')
pylab.legend(loc=0)
