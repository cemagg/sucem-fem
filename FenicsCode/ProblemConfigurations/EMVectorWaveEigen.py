import dolfin
from FenicsCode import Forms 
from FenicsCode import Materials 
from FenicsCode import BoundaryConditions 
from FenicsCode import SystemMatrices

from scipy.sparse.linalg.eigen.arpack import speigs
import FenicsCode.BoundaryConditions.essential
import FenicsCode.BoundaryConditions.container

class CombineForms(Forms.CombineGalerkinForms):
    def get_forms(self):
        m = self.interior_forms.get_mass_form()
        s = self.interior_forms.get_stiffness_form()
        # I suspect that using e.g. the bilinear part of a first order
        # ABC in an eigen solution could potentially be used for eigen
        # solution, so I'm including the bilinear boundary condition
        # forms here
        BC_bilin_forms = self.boundary_conditions.get_bilinear_form()
        return dict(M=m, S=s, BC=BC_bilin_forms)

class EigenProblem(object):
    element_type = "Nedelec 1st kind H(curl)"
    CombineForms = CombineForms        
    material_regions = None
    region_meshfunction = None
    bcs = None
    
    def set_mesh(self, mesh):
        self.mesh = mesh

    def set_basis_order(self, order):
        self.basis_order = order

    def set_material_regions(self, material_regions):
        """Set material region properties

        See documentation of Materials.MaterialPropertiesFactory for input format
        """
        self.material_regions = material_regions

    def set_region_meshfunction(self, region_meshfunction):
        self.region_meshfunction = region_meshfunction

    def init_problem(self):
        self.function_space = dolfin.FunctionSpace(
            self.mesh, self.element_type, self.basis_order)
        self.material_properties = self._get_material_properties()
        self.material_functions = self._get_material_functions()
        self.interior_forms = self._get_interior_forms()
        self.boundary_conditions = self._get_boundary_conditions()
        self.combined_forms = self._get_combined_forms()
        self.system_matrices = self._get_system_matrices()
        
    def _get_material_properties(self):
        mat_props_fac = Materials.MaterialPropertiesFactory(
            self.material_regions)
        return mat_props_fac.get_material_properties()

    def _get_material_functions(self):
        mat_func_fac = Materials.MaterialFunctionFactory(
            self.material_properties, self.region_meshfunction, self.mesh)
        return mat_func_fac.get_material_functions('eps_r', 'mu_r')

    def _get_interior_forms(self):
        int_form = Forms.EMGalerkinInteriorForms()
        int_form.set_material_functions(self.material_functions)
        int_form.set_function_space(self.function_space)
        return int_form

    def _get_boundary_conditions(self):
        if self.bcs is None:
            self.bcs = FenicsCode.BoundaryConditions.container.BoundaryConditions()
            bc = FenicsCode.BoundaryConditions.essential.EssentialBoundaryCondition()
            bc.set_function_space(self.function_space)
            class bcSubDomain(dolfin.SubDomain):
                def inside(self, x, on_boundary):
                    return on_boundary
            bc.init_with_subdomain(bcSubDomain(), 0)
            self.bcs.add_boundary_condition(bc)
        return self.bcs

    def _get_combined_forms(self):
        comb_forms = self.CombineForms()
        comb_forms.set_interior_forms(self.interior_forms)
        comb_forms.set_boundary_conditions(self.boundary_conditions)
        return comb_forms
            
    def _get_system_matrices(self):
        bilin_forms = self.combined_forms.get_forms()
        sysmats = SystemMatrices.SystemMatrices()
        sysmats.set_matrix_forms(bilin_forms)
        sysmats.set_boundary_conditions(self.boundary_conditions)
        return sysmats.calc_system_matrices()

class DefaultEigenSolver(object):
    def set_eigenproblem(self, eigenproblem):
        """Sets initialised instance of EigenProblem to solve"""
        self.eigenproblem = eigenproblem

    def set_sigma(self, sigma):
        """Spectrum shift (sigma) to apply to k^2"""
        self.sigma = sigma

    def solve_problem(self, nev, ncv=None):
        """Solve problem for eigenvalues and eigenvectors
        Input Values
        -------------
        @param nev: Number of eigenpairs to compute
        @param ncv: Number of Arnoldi basisvectors to use. If None, default to 2*nev+1
        Return Values
        -------------
        (eigs_w, eigs_v) with

        eigs_w -- array of n_eig eigen values corresponding to k^2 for the
            mode's resonant wavenumber

        eigs_v -- 2D array ofshape n_eig x problem_dim eigen vectors
            corresponding to the modal distributions of the eigenvalues.
            eigs_v[i] is the eigenvector corresponding to eigs_w[i]
        
        """
        M = self.eigenproblem.system_matrices['M']
        S = self.eigenproblem.system_matrices['S']
        solve_mat = S - self.sigma*M
        lu = dolfin.LUSolver(solve_mat)
        lu.parameters["reuse_factorization"] = True
        lu.parameters["report"] = False
        bb = dolfin.Vector(M.size(0))
        xx = dolfin.Vector(M.size(0))
        def sigma_solve(b):
            bb[:] = b
            lu.solve(xx, bb)
            return xx[:]
        M_matvec = lambda x: M*x
        eigs_w, eigs_v = speigs.ARPACK_gen_eigs(
            M_matvec, sigma_solve, M.size(0), self.sigma, nev, ncv=ncv)

        return eigs_w, eigs_v.T
