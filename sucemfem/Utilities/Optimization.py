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

import dolfin

def set_dolfin_optimisation(enable=True):
    """Set dolfin optimization options
    
    @keyword enable:  
        True -> optimization enabled
        False -> disabled
        (default: True).
    """
    dolfin.parameters['optimize_form'] = enable
    dolfin.parameters['optimize'] = enable
    dolfin.parameters['optimize_use_dofmap_cache'] = enable
    dolfin.parameters['optimize_use_tensor_cache'] = enable
    dolfin.parameters['form_compiler']['optimize'] = enable
    dolfin.parameters['form_compiler']['cpp_optimize'] = enable

