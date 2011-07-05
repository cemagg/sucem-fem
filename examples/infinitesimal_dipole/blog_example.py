# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import numpy as N
import dolfin 
from dolfin import dot, cross, curl, inner, dx, ds
import scipy.sparse

# Some constants
pi = N.pi
e = N.e
# See http://en.wikipedia.org/wiki/Speed_of_light
mu0 = N.pi*4*1e-7                       # Per definition 
c0 = 299792458                          # Per definition
eps0 = 10**7/4/N.pi/c0/c0               # Per definition
# Per definition, see http://en.wikipedia.org/wiki/Impedance_of_free_space#Relation_to_other_constants
Z0 = mu0*c0                              

# Define problem parameters
eps_r = 1;
mu_r = 1;
freq = 1e9;
lam = c0/freq
source_coord = N.array([0,0,0], dtype=N.float64)
source_point = dolfin.Point(*source_coord)
# Field evaluation points, from 1/10 of a wavelength to 0.99 of a
# wavelength from the dipole
field_pts = N.array([lam,0,0])*(N.arange(90)/100+1/10)[:, N.newaxis]
I = 1                                   # Dipole current
l = lam/1000                            # Dipole length
source_value = N.array([0,0,1.])*I*l

# Solution parameters
order = 2                               # Basis element order
max_edge_len = lam/6
domain_size = N.array([2*lam]*3)
solver = 'iterative'                    # Solver type, iterative or direct

def dolfin_ublassparse_to_scipy_csr ( A, dtype=None, imagify=False ):
    """
    convert a DOLFIN uBLASSparseMatrix to a scipy.sparse.csr_matrix()
    
    @param A: a DOLFIN uBLASSparseMatrix
    @param dtype: the numpy data type to use to store the matrix
    @param imagify: multiply the original matrix data by 1j
    """
    import scipy.sparse
    # get the sparse data from the input matrix
    (row,col,data) = A.data()   # get sparse data
    col = N.intc(col)
    row = N.intc(row)
    n = A.size(0)
    if imagify: data = data*1j
    A_sp = scipy.sparse.csr_matrix( (data,col,row), shape=(n,n), dtype=dtype)
    
    return A_sp

def solve_sparse_system ( A, b ):
    """
    This function solves the sparse linear system Ax = b for A a scipy sparse matrix, and b a numpy or scipy array
    
    A inverse diagonal preconditioner is used with the bicgstab iterative solver
    
    @param A: a square matrix 
    @param b: the RHS vector
    """
    import scipy.sparse
    import scipy.sparse.linalg
    
    x, info = scipy.sparse.linalg.bicgstab(
        A, b, M=scipy.sparse.spdiags(1./A.diagonal(), 0, A.shape[0], A.shape[1]) )
    
    assert ( info == 0 )
    
    return x

def calc_pointsource_contrib(V, source_coords, source_value):
    """Calculate the RHS contribution of a current point source (i.e. electric dipole)
    Input Values
    -------------
    @param V: dolfin FunctionSpace object
    @param source_coords: length 3 array with x,y,z coordinates of point source
    @param source_value: length 3 array with x,y,z componets of source current 
    Return Values
    -------------
    (dofnos, rhs_contribs) with

    dofnos -- Array of degree of freedom indices of the source contribution

    rhs_contribs -- Numerical values of RHS contribution, such that
        RHS[dofnos] += rhs_contribs will add the current source to the system.

    """
    source_coords = N.asarray(source_coords, dtype=N.float64)
    source_value = N.asarray(source_value)
    dm = V.dofmap()
    dofnos = N.zeros(dm.max_cell_dimension(), dtype=N.uint32)
    source_pt = dolfin.Point(*source_coords)
    try:
        cell_index = V.mesh().any_intersected_entity(source_pt)
    except StandardError:
        # CGAL as used by dolfin to implement intersection searches
        # seems to break with 1-element meshes
        if dolfin.Cell(V.mesh(), 0).intersects(source_pt):
            cell_index = 0
        else: raise
    c = dolfin.Cell(V.mesh(), cell_index)
    # Check that the source point is in this element    
    assert(c.intersects_exactly(source_pt)) 

    dm.tabulate_dofs(dofnos,  c)
    finite_element = V.dolfin_element()
    no_basis_fns = finite_element.space_dimension()
    # Vector valued elements have rank of 1
    assert(finite_element.value_rank() == 1)
    # Vector valued elements have only one rank (i.e. 0) along which
    # dimensions are defined. This is the dimension that the basis
    # function value vector is. Since we have 3D Nedelec elements here
    # this should be 3
    bf_value_dimension = finite_element.value_dimension(0)
    el_basis_vals = N.zeros((no_basis_fns, bf_value_dimension), dtype=N.float64)
    finite_element.evaluate_basis_all(el_basis_vals, source_coords, c)
    # Compute dot product of basis function values and source value
    rhs_contribs = N.sum(el_basis_vals*source_value, axis=1)
    return dofnos, rhs_contribs

def main():

    # Define mesh
    domain_subdivisions = N.array(N.ceil(N.sqrt(2)*domain_size/max_edge_len), N.uint)
    print 'Numer of domain subdomain_subdivisions: ', domain_subdivisions

    mesh = dolfin.UnitCube(*domain_subdivisions)
    # Transform mesh to correct dimensions
    mesh.coordinates()[:] *= domain_size
    # Centred around [0,0,0]
    mesh.coordinates()[:] -= domain_size/2
    ## Translate mesh slightly so that source coordinate lies at
    ## centroid of an element
    source_elnos = mesh.all_intersected_entities(source_point)
    closest_elno = source_elnos[(N.argmin([source_point.distance(dolfin.Cell(mesh, i).midpoint())
                                  for i in source_elnos]))]
    centre_pt = dolfin.Cell(mesh, closest_elno).midpoint()
    centre_coord = N.array([centre_pt.x(), centre_pt.y(), centre_pt.z()])
    # There seems to be an issue with the intersect operator if the
    # mesh coordinates are changed after calling it for the first
    # time. Since we called it to find the centroid, we should init a
    # new mesh
    mesh_coords = mesh.coordinates().copy()
    mesh = dolfin.UnitCube(*domain_subdivisions)
    mesh.coordinates()[:] = mesh_coords
    mesh.coordinates()[:] -= centre_coord
    ##

    # Define function space
    V = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", order)

    k_0 = 2*N.pi*freq/c0                    # Freespace wave number
    # Definite test- and trial functions
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)

    # Define the bilinear forms
    m = eps_r*inner(v, u)*dx                # Mass form
    s = (1/mu_r)*dot(curl(v), curl(u))*dx   # Stiffness form

    n = V.cell().n                           # Get the surface normal
    s_0 = inner(cross(n, v), cross(n, u))*ds # ABC boundary condition form

    # Assemble forms using uBLASS matrices so that we can easily export to scipy
    print 'Assembling forms'
    M = dolfin.uBLASSparseMatrix()
    S = dolfin.uBLASSparseMatrix()
    S_0 = dolfin.uBLASSparseMatrix()
    dolfin.assemble(m, tensor=M, mesh=mesh)
    dolfin.assemble(s, tensor=S, mesh=mesh)
    dolfin.assemble(s_0, tensor=S_0, mesh=mesh)
    print 'Number of degrees of freedom: ', M.size(0)


    # Set up RHS
    b = N.zeros(M.size(0), dtype=N.complex128)
    dofnos, rhs_contrib = calc_pointsource_contrib(V, source_coord, source_value)

    rhs_contrib = 1j*k_0*Z0*rhs_contrib
    b[dofnos] += rhs_contrib

    Msp = dolfin_ublassparse_to_scipy_csr(M)
    Ssp = dolfin_ublassparse_to_scipy_csr(S)
    S_0sp = dolfin_ublassparse_to_scipy_csr(S_0)

    # A is the system matrix that must be solved
    A = Ssp - k_0**2*Msp + 1j*k_0*S_0sp     


    solved = False;
    if solver == 'iterative':
        # solve using scipy bicgstab
        print 'solve using scipy bicgstab'
        x = solve_sparse_system ( A, b )
    elif solver == 'direct':
        import scipy.sparse.linalg
        A_lu = scipy.sparse.linalg.factorized(A.T)
        x = A_lu(b)
    else: raise ValueError("solver must have value 'iterative' or 'direct'")

    dolfin.set_log_active(False) # evaluation seems to make a lot of noise
    u_re = dolfin.Function(V)
    u_im = dolfin.Function(V)
    # N.require is important, since dolfin seems to expect a contiguous array
    u_re.vector()[:] = N.require(N.real(x), requirements='C')
    u_im.vector()[:] = N.require(N.imag(x), requirements='C')
    E_field = N.zeros((len(field_pts), 3), dtype=N.complex128)
    for i, fp in enumerate(field_pts):
        try: E_field[i,:] = u_re(fp) + 1j*u_im(fp)
        except (RuntimeError, StandardError): E_field[i,:] = N.nan + 1j*N.nan

    import pylab as P
    r1 = field_pts[:]/lam
    x1 = r1[:,0]
    E_ana = N.abs(analytical_result)
    E_num = E_field
    P.figure()
    P.plot(x1, N.abs(E_num[:,0]), '-g', label='x_num')
    P.plot(x1, N.abs(E_num[:,1]), '-b', label='y_num')
    P.plot(x1, N.abs(E_num[:,2]), '-r', label='z_num')
    P.plot(analytical_pts, E_ana, '--r', label='z_ana')
    P.ylabel('E-field Magnitude')
    P.xlabel('Distance (wavelengths)')
    P.legend(loc='best')
    P.grid(True)
    P.show()



## Analytical result z-component. x- and y- components are zero for
## z-directed dipole.
analytical_result = N.array([- 2.428440766843998, - 2.386808396555343, - 2.341640195751903, 
- 2.2930497650128, - 2.24115911023354, - 2.186098288024988, 
- 2.128005028098966, - 2.067024333820723, - 2.003308062175806, 
- 1.937014484462164, - 1.868307829077282, - 1.797357807824657, 
- 1.724339127213784, - 1.649430986272914, - 1.572816562433923, 
- 1.494682487083704, - 1.415218312406344, - 1.33461597116493, 
- 1.253069231091067, - 1.170773145564037, - 1.087923502269868, 
- 1.004716271533579, - 0.92134705601524, - 0.83801054345257, 
- 0.75489996411929, - 0.67220655464987, - 0.59011902985691, - 0.5088230641385, 
- 0.42850078403837, - 0.34933027348209, - 0.27148509316891, 
- 0.19513381554938, - 0.12043957676596, - 0.047559646875914, 0.02335498038611, 
0.092159977113948, 0.15871804888206, 0.2228992794945, 0.28458145185751, 
0.34365034408151, 0.4, 0.45353297337717, 0.50416054516209, 0.55180291323563, 
0.59638935418576, 0.63785835673786, 0.67615772655851, 0.71124466224343, 
0.74308580239326, 0.77165724377357, 0.79694453064779, 0.81894261546378, 
0.83765579116569, 0.8530975954925, 0.86529068771273, 0.87426669833161, 
0.88006605239108, 0.88273776706514, 0.88233922433266, 0.87893591958607, 
0.87260118710818, 0.86341590341905, 0.85146816956193, 0.83685297345958, 
0.8196718335315, 0.80003242481698, 0.77804818889933, 0.7538379289727, 
0.72752539143374, 0.69923883541729, 0.66911059172678, 0.63727661263676, 
0.60387601406694, 0.56905061164349, 0.53294445217542, 0.49570334208004, 
0.45747437429304, 0.4184054551951, 0.37864483307813, 0.33834062966058, 
0.29764037614242, 0.25669055526701, 0.21563615082877, 0.17462020603233, 
0.1337833920718, 0.093263588256616, 0.053195474964451, 0.013710140651745, 
- 0.025065295901765, - 0.063008046995554]) + 1j*N.array(
[13.67055420062157, 10.08025338822953, 7.650988748140257, 
5.963849498300017, 4.767565345646867, 3.905342078904754, 3.275987685057256, 
2.812270066294151, 2.468375209936027, 2.212381731606761, 2.021605548302317, 
1.879641375967055, 1.774436926405598, 1.697012096444633, 1.640590472571951, 
1.6, 1.571252742280137, 1.551245884852265, 1.537546133017615, 
1.528232311858383, 1.521779134146906, 1.516970450175031, 1.512833853790133, 
1.50859092358448, 1.503619024164276, 1.49742173323039, 1.489605760148637, 
1.479862788775132, 1.467955083455426, 1.453703990823293, 1.436980684355172, 
1.417698656388281, 1.395807579364665, 1.371288245588828, 1.344148360710643, 
1.314419016143676, 1.282151703799558, 1.247415765842435, 1.210296194829229, 
1.170891717214529, 1.129313106962401, 1.085681686818093, 1.040127983325753, 
0.9927905084448, 0.94381464600689, 0.89335162556531, 0.84155756964876, 
0.78859260321841, 0.73462001638011, 0.67980547322871, 0.62431626118421, 
0.5683205763861, 0.51198684169512, 0.45548305465299, 0.39897616340345, 
0.34263146910946, 0.28661205383211, 0.23107823318534, 0.17618703336003, 
0.12209169233316, 0.068941185252098, 0.016879774117985, - 0.033953418007558, 
- 0.083424807977039, - 0.13140672837035, - 0.17777777777778, 
- 0.22242314517492, - 0.26523490801796, - 0.3061123037067, - 0.34496197409247, 
- 0.38169818274731, - 0.41624300475771, - 0.4485264888595, - 0.47848679178995, 
- 0.50607028479605, - 0.53123163230589, - 0.55393384283944, 
- 0.57414829230784, - 0.59185471992389, - 0.60704119702119, - 0.6197040691548, 
- 0.62984787193106, - 0.63748522108931, - 0.64263667743131, 
- 0.64533058726692, - 0.64560289911492, - 0.64349695746642, 
- 0.63906327448423, - 0.63235928057537, - 0.6234490548339])
# x-coordinate of analytical solution
analytical_pts = N.array([0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 
0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 
0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 
0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 
0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 
0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 
0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]) 

if __name__ == '__main__':
    main()
