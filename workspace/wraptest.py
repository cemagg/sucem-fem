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
# Authors
# Neilen Marais <nmarais@gmail.com>

import os
from scipy import *
from femfeko import *

def get_dofmap(map):
    n = coupled_td.get_n(map)
    return coupled_td.wrap_dofmap(n, map)
    
femfeko.run_femfeko()

# n = coupled_td.get_n('e')
# print coupled_td.wrap_dofmap(n, 'e')
# n = coupled_td.get_n('d')
# print coupled_td.wrap_dofmap(n, 'd')
# n = coupled_td.get_n('h')
# print coupled_td.wrap_dofmap(n, 'h')
# n = coupled_td.get_n('b')
# print coupled_td.wrap_dofmap(n, 'b')

geomwrap.init_geom()

print geomwrap.vertex_coords
