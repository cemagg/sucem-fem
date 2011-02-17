import numpy as N
from scipy import sparse, cross
from itertools import izip
from warnings import warn

from NewCode.Exceptions import AllZero

numint_tol_fac = 111
eps = N.finfo(float).eps*numint_tol_fac

def local_self_projection_matrix(el):
    physVals = el.physVals()
    intg = el.rule.integrateFun
    # Calculate the dot product of two functions at each integration point
    # using the sum() call, then integrate to build the local matrix
    local_mat = N.array([[intg(N.sum(fn_i*fn_j, axis=1))
                          for fn_j in physVals]
                         for fn_i in physVals], N.float64)
    local_mat[N.abs(local_mat)/N.max(N.abs(local_mat)) < eps] = 0.
    local_mat *= el.size          # Integration assumes unit size 
    return local_mat
    

def local_boundary_matrix(el, faces):
    boundarySet = faces.boundarySet
    no_dofs = el.noDOFs.element
    local_mat = N.zeros((no_dofs, no_dofs), N.float64)
    intg = el.faceRule.integrateFun
    face_intg_pts = el.faceRule.evalPoints()
    normals = el.outwardNormals()
    
    for l_faceno, faceno in enumerate(el.facenos):
        if faceno not in boundarySet: continue
        normal = normals[l_faceno]
        face = faces[faceno]
        area = face.area()
        tet_intg_pts = N.array(
            [el.face_coords2vol_coords(l_faceno, pt) for pt in face_intg_pts])
        physVals = el.physValsAtPoints(tet_intg_pts)
        for bf_i in range(len(physVals)):
            physVals[bf_i,:] = [cross(normal, val) for val in physVals[bf_i]]
        for i, fn_i in enumerate(physVals):
            for j, fn_j in enumerate(physVals):
                local_mat[i,j] += area*intg(N.sum(fn_i*fn_j, axis=1))
        
    return local_mat

class ConstructionMatrix(object):
    def __init__(self, shape=None, dtype=N.float64):
        self.dtype = dtype ; self.shape = shape
        self.row = [] ; self.col = []
        self.data = [] 

    def insert(self, locmat, perm_i, perm_j=None):
        """Construct global row/col/data given local matrices and permutations

        Data can be combined into a COO matrix. Repeated matrix indices will be added
        """
        if perm_j is None: perm_j = perm_i
        row, col, data = self.row, self.col, self.data
        for l_i, g_i in izip(*perm_i):
            for l_j, g_j in izip(*perm_j):
                d = locmat[l_i,l_j]
                if d != 0:
                    row.append(g_i)
                    col.append(g_j)
                    data.append(d)

    def get_matrix(self):
        return sparse.coo_matrix((self.data, (self.row, self.col)),
                                 dtype=self.dtype, shape=self.shape)

class DiagonalConstructionMatrix(object):
    def __init__(self, shape=None, dtype=N.float64):
        assert len(shape) == 2
        assert shape[0] == shape[1]
        size = shape[0]
        self.dtype = dtype ; self.shape = shape
        self.data = N.zeros(size, dtype=dtype)

    def insert(self, locmat, perm_i):
        loc, glob = perm_i
        self.data[glob] += locmat[loc, loc]

    def get_matrix(self):
        return sparse.dia_matrix(
            ([self.data], 0), shape=self.shape, dtype=self.dtype)

def self_projection_matrix(disc, dtype=N.float64):
    """
    Calculate the mass matrix of system disc, returning a scipy.sparse.lil_matrix
    """
    no_dofs = disc.totalDOFs
    diag = disc.diagonalMetric
    CM = DiagonalConstructionMatrix if diag else ConstructionMatrix
    proj_constr = CM(shape=(no_dofs, no_dofs), dtype=dtype)

    lm_uncalced = True ; different = not disc.identicalElements
    for el in disc.elements:
        try: perm = el.permutation()
        except AllZero: continue
        if lm_uncalced or different: local_mat = local_self_projection_matrix(el)
        lm_uncalced = False
        proj_constr.insert(local_mat, el.permutation())
    return proj_constr.get_matrix()

def boundary_matrix(disc, dtype=N.float64):
    """
    Calculate the boundary surface matrix of system disc, returning a scipy.sparse.lil_matrix
    """
    no_dofs = disc.totalDOFs
    bdry_mat = sparse.lil_matrix(shape=(no_dofs, no_dofs), dtype=dtype)
    for el in disc.elements:
        local_mat = local_boundary_matrix(el, disc.mesh.faces)
        insert_global(bdry_mat, local_mat, el.permutation())
    return bdry_mat

def local_projection_matrix(elA, elB):
    physValsA = elA.physVals()
    physValsB = elB.physVals()
    intg = elA.rule.integrateFun
    # The integration orders must be the same for both elements
    assert(elA.rule.order == elB.rule.order) 
    # Calculate the dot product of two functions at each integration point
    # using the sum() call, then integrate to build the local matrix
    local_mat = N.array([[intg(N.sum(fn_i*fn_j, axis=1))
                              for fn_j in physValsB]
                         for fn_i in physValsA], N.float64)
    local_mat[N.abs(local_mat)/N.max(N.abs(local_mat)) < eps] = 0.
    local_mat *= elA.size          # Integration assumes unit size
    return local_mat

def projection_matrix(discA, discB=None, dtype=N.float64):
    """
    Calculate the projection matrix of two discretisers, or self-projection for one discretiser.

    Returning a scipy.sparse.lil_matrix
    """
    if discB is None: return self_projection_matrix(discA, dtype=dtype)
    if hasattr(discA.matrix, 'compoundProjection'):
        return discA.matrix.compoundProjection(discB).T
    no_dofs_i, no_dofs_j = discA.totalDOFs, discB.totalDOFs
    assert(discA.mesh is discB.mesh)
    proj_constr = ConstructionMatrix(shape=(no_dofs_i, no_dofs_j), dtype=dtype)
    lm_uncalced = True
    different = not (discA.identicalElements and discB.identicalElements)
    for (elA, elB) in izip(discA.elements, discB.elements):
        try: permA, permB = elA.permutation(), elB.permutation()
        except AllZero:     
            continue                    # If either element has all zero values.
        if lm_uncalced or different: local_mat = local_projection_matrix(elA, elB)
        lm_uncalced = False
        proj_constr.insert(local_mat, permA, permB)
    return proj_constr.get_matrix()


def lumpy_projection_matrix(discA, discB, dtype=N.float64):
    """
    discA is the discretiser being projected onto, discB is the \"master\"
    discretiser that defines the integration rule
    """
    no_dofs_i, no_dofs_j = discA.totalDOFs, discB.totalDOFs
    assert(discA.mesh is discB.mesh)
    proj_constr = ConstructionMatrix(shape=(no_dofs_i, no_dofs_j), dtype=dtype)
    coord = {(1.,0.,0.):0,
             (0.,1.,0.):1,
             (0.,0.,1.):2}
    for (elA, elB) in izip(discA.elements, discB.elements):
        try: permA, permB = elA.permutation(), elB.permutation()
        except AllZero:     
            continue                    # If either element has all zero values.
        physValsA = elA.alt_physVals()
        physValsB = elB.physVals()
        intg = elB.rule.integrateFun
        local_mat = N.array([[intg(
            N.sum(fn_i[coord[tuple(bf.direction)]]*fn_j, axis=1))
                              for fn_j, bf in izip(physValsB, elB.basisSet)]
                             for fn_i in physValsA], N.float64)
        local_mat[N.abs(local_mat)/N.max(N.abs(local_mat)) < eps] = 0.
        local_mat *= elA.size          # Integration assumes unit size
        proj_constr.insert(local_mat, permA, permB)
    return proj_constr.get_matrix()
    
    
def insert_global(glob_mat, local_mat, perm_i, perm_j=None):
    """Add local_mat contributions to glob_mat"""
    warn("Use ConstructionMatrix method for MUCH faster global matrix construction")
    if perm_j is None: perm_j = perm_i
    for local_i, glob_i in izip(*perm_i):
        for local_j, glob_j in izip(*perm_j):
            glob_mat[glob_i, glob_j] += local_mat[local_i, local_j]
            
def set_global(glob_mat, local_mat, perm_i, perm_j=None):
    """Set glob_mat entries equal to local_mat entries

    Differs from insert_global in that the global matrix entry is replaced by
    the local matrix entry rather than having the local matrix entry added.
    I.e.

    glob_mat[glob_i, glob_j] = local_mat[local_i, local_j]
    
    """
    if perm_j is None: perm_j = perm_i
    for local_i, glob_i in izip(*perm_i):
        for local_j, glob_j in izip(*perm_j):
            glob_mat[glob_i, glob_j] = local_mat[local_i, local_j]
            
