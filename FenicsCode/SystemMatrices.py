from __future__ import division

import dolfin 
from FenicsCode import Forms

class SystemMatrices(object):
    MatrixClass = dolfin.PETScMatrix

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
