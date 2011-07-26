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
import numpy as N

pi = N.pi
e = N.e
# See http://en.wikipedia.org/wiki/Speed_of_light
mu0 = N.pi*4*1e-7                       # Per definition 
c0 = 299792458                          # Per definition
eps0 = 10**7/4/N.pi/c0/c0               # Per definition
# Per definition, see http://en.wikipedia.org/wiki/Impedance_of_free_space#Relation_to_other_constants
Z0 = mu0*c0                              

c = mu = eps = 1.                       # Our simplified choice

# stability timestep for lumped unit cuboid hexahedra to order. Actual
# stability will be given by dt = h*stability_factor, where h is the gridstep
# size.
lumped_stability_factors = {
    1:1/N.sqrt(3),
    2:2/N.sqrt(72),
    3:2/N.sqrt(222.9329665284213),
    4:2/N.sqrt(550.04544985070766),
    5:2/N.sqrt(1175.8767213081619),
    6:2/N.sqrt(2248.6728409490497),
    7:2/N.sqrt(73943.67705255)}
