# Authors:
# Neilen Marais <nmarais@gmail.com>
import dolfin as dol
from dolfin import dot, curl, inner, dx
from scipy.sparse.linalg.eigen.arpack import speigs

# Test for PETSc and SLEPc
if not dol.has_la_backend("PETSc"):
    print "DOLFIN has not been configured with PETSc. Exiting."
    exit()

if not dol.has_slepc():
    print "DOLFIN has not been configured with SLEPc. Exiting."
    exit()

class CavityDims(object): pass
cdims = CavityDims()
cdims.a, cdims.b, cdims.c = 29,23,19


# Define mesh, function space
mesh = dol.UnitCube(1,1,1)
mesh.coordinates()[:] *= [cdims.a,cdims.b,cdims.c]
V = dol.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", 4)

# Define basis and bilinear form
u = dol.TrialFunction(V)
v = dol.TestFunction(V)
m = inner(v, u)*dx                      # Mass form
s = dot(curl(v), curl(u))*dx            # Stiffness form

# Assemble smass form
M = dol.PETScMatrix()
S = dol.PETScMatrix()
dol.assemble(m, tensor=M, mesh=mesh)
dol.assemble(s, tensor=S, mesh=mesh)

sigma = 0.03
smat = S - sigma*M
#lu = dol.LUSolver(S - sigma*M)
lu = dol.LUSolver(smat)
lu.parameters["reuse_factorization"] = True
lu.parameters["report"] = False
bb = dol.Vector(M.size(0))
xx = dol.Vector(M.size(0))
def sigma_solve(b):
    bb[:] = b
    lu.solve(xx, bb)
    return xx[:]
M_matvec = lambda x: M*x

arpack_eigs,v = speigs.ARPACK_gen_eigs(M_matvec, sigma_solve, M.size(0), sigma, 51, ncv=91)

# Create eigensolver
esolver = dol.SLEPcEigenSolver(S,M)
esolver.parameters["spectrum"] = "smallest real"
esolver.parameters["solver"] = "arnoldi"
esolver.parameters["spectral_shift"] = sigma
esolver.parameters["spectral_transform"] = "shift-and-invert"
esolver.solve()

eigs = [esolver.get_eigenvalue(i)[0] for i in range(4)]
filtered_eigs = []
for i in range(S.size(0)):
    ev_r, ev_i = esolver.get_eigenvalue(i)
    if ev_r > 0.001:
        filtered_eigs.append(ev_r)


print sorted(arpack_eigs)[0:4]
print eigs[0:4]
print filtered_eigs[0:4]

    


