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
# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import numpy as N
from postproc_code.AnalyticResults import PEC_cavity, err_percentage
from NewCode import Utilities

class EigenErrors(object):
    pass

def calc_errs(eigenvalues, geomname='rect-white', RMS_modes=4):
    errs = EigenErrors()
    errs.eigenvalues = eigenvalues
    errs.geomname = geomname
    errs.RMS_modes = RMS_modes 
    errs.ares = PEC_cavity[geomname]
    errs.eigenvalue_errors = err_percentage(errs.ares, eigenvalues)
    errs.RMS_eigenvalue_errors = Utilities.RMS(errs.eigenvalue_errors[0:4])
    return errs

def print_errs(errs):
    print 'eigenvalues: \n', errs.eigenvalues
    print 'RMS error (%%) over first %d modes: %s' % (
        errs.RMS_modes, str(errs.RMS_eigenvalue_errors))
    print 'mode errors (%): \n', errs.eigenvalue_errors
