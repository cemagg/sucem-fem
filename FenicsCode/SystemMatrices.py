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

import dolfin 
from FenicsCode import Forms

class SystemMatrices(object):
    MatrixClass = dolfin.PETScMatrix

    def set_matrix_class(self, matrix_class):
        """Set matrix class to use for system matrix.

        Overrides the default value (dolfin.PETScMatrix) set in the
        class definition
        
        """
        self.MatrixClass = matrix_class

    def set_matrix_forms(self, matrix_forms):
        """Set matrix_forms with a dict mapping matrix names to bilinear forms"""
        self.matrix_forms = matrix_forms

    def set_boundary_conditions(self, boundary_conditions):
        """Set boundary_conditions with instance of BoundaryConditions"""
        self.boundary_conditions = boundary_conditions

    def calc_system_matrices(self):
        """Calculate and return system matrices in a dict"""
        system_matrices = dict()
        for matname, form in self.matrix_forms.items():
            mat = self.MatrixClass()
            if isinstance(form, Forms.NullForm):
                mat = None
            else:
                dolfin.assemble(form, tensor=mat)
                self.boundary_conditions.apply_essential(mat)
            system_matrices[matname] = mat

        return system_matrices
            
            
class SystemVectors(object):
    VectorClass = dolfin.Vector
    def set_vector_forms(self, vector_forms):
        """Set vector_forms with a dict mapping vector names to linear forms"""
        self.vector_forms = vector_forms

    def set_boundary_conditions(self, boundary_conditions):
        """Set boundary_conditions with instance of BoundaryConditions"""
        self.boundary_conditions = boundary_conditions

    def calc_system_vectors(self):
        """Calculate and return system vectors in a dict"""
        system_vectors = dict()
        for vecname, form in self.vector_forms:
            vec = self.VectorClass()
            dolfin.assemble(form, tensor=vec)
            self.boundary_conditions.apply_essential(vec)
            system_vectors[vecname] = vec

        return system_vectors
