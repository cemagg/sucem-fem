# Authors:
# Evan Lezar <mail@evanlezar.com>
"""
    Work on implementing the Pelosi formulation of an h-plane waveguide
"""
import numpy as N
import pylab as P
from dolfin import *
import scipy.sparse
import csv

import sys
sys.path.insert(0, '../../')
from FenicsCode.Utilities.LinalgSolvers import solve_sparse_system
from FenicsCode.Utilities.Converters import dolfin_ublassparse_to_scipy_csr
del sys.path[0]
#from solvers import CutoffSquared
#parameters["linear algebra backend"] = "uBLAS"

a = 18.35
b = 9.175
l = 10.0

Z_0 = 377.0

def plot_two_axis_data(X, Y1_data, Y2_data, Y_labels = ['y1', 'y2'], Y1_styles = None, Y2_styles = None, Y1_lables = None, Y2_labels = None, Y1_limits = None, Y2_limits = None):
    
    def plot_axis_data(ax, X, data_set, data_styles, data_labels, ylimit):
        for d in range(len(data_set)):
            data = data_set[d]
            try:
                style = data_styles[d]
            except:
                style = ''
            
            try:
                labels = data_labels[d]
            except:
                labels = None
                
            ax.plot(X, data, style, label=labels, linewidth=2)
            
            if not ylimit == None:
                P.ylim(ylimit)
    
        
    fig = P.figure()
    ax1 = fig.add_subplot(111)
    plot_axis_data(ax1, X, Y1_data, Y1_styles, Y1_lables, Y1_limits)
    P.ylabel(Y_labels[0])
    P.xticks([])
    ax2 = ax1.twinx()
    plot_axis_data(ax2, X, Y2_data, Y2_styles, Y2_labels, Y2_limits)
    P.ylabel(Y_labels[1])


def phase(a, adjust_range = True, perform_mod = True):
    results = N.arctan2(a.imag, a.real)*180/N.pi 
    
    if adjust_range:
        if N.isscalar(results):
            if results < 0:
                results = 360 - results
            if perform_mod:
                return N.mod(results, 360)
            else:
                return results
        
        for i in range(len(results)):
            if results[i] < 0:
                results[i] = 360 + results[i]
        if perform_mod: 
            results = N.mod(results, 360)
        else:
            return results
    
    return results 

def magnitude(a, dB=True):
    tmp = N.abs(a)
    if dB:
        return 20*N.log10(tmp)
    else:
        return tmp
    
def k_o(f):
    return 2*N.pi*f/3e8 

def Beta_km(f, m = 1):
    k_tm = m*N.pi/a
    k_o_local = k_o(f)

    if k_o_local*k_o_local >= k_tm*k_tm:
        return N.sqrt(k_o_local*k_o_local - k_tm*k_tm)
    else:
        return -1j*N.sqrt(k_tm*k_tm - k_o_local*k_o_local)


class WaveguideWalls(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary:
            return (x[1] <= DOLFIN_EPS) or (x[1] >= (a - DOLFIN_EPS))
        return False

#class WaveguideWalls(SubDomain):
#    def inside(self, x, on_boundary):
#        c = 4.587
#        w = 9.175
#        d = 1.0
#        r = False
#        if ((x[1] <= DOLFIN_EPS) or (x[1] >= ( a -DOLFIN_EPS)))  and ( x[0] > DOLFIN_EPS) and ( x[0] < ( l - DOLFIN_EPS) ):
#            r = True
#        elif ((x[0] >= d - DOLFIN_EPS ) and (x[0] <= d + 0.5 + DOLFIN_EPS)):
#            if (x[1] <= (c + DOLFIN_EPS)) or (x[1] >= (c + w - DOLFIN_EPS)): 
#                r = True
#            else:
#                r = False 
#        else:
#            r = False
#    
#        if r:
#            P.plot([x[0],], [x[1],], 'r.')
#        else:
#            P.plot([x[0],], [x[1],], 'b.')
#            
#        return r

#class WaveguideWalls(SubDomain):
#    def inside(self, x, on_boundary):
#        if on_boundary:
#            return (x[0] > DOLFIN_EPS) and (x[0] < (l))
#        return False

class InputPort(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] <= DOLFIN_EPS

class OutputPort(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] >= (l -DOLFIN_EPS)


def reference_phase(f_range, adjust_range, perform_mod):
    bb = N.zeros(len(f_range), dtype=N.complex128)

    for i in range(len(f_range)):
        bb[i] = Beta_km(f_range[i]*f_scale)
        
    #print bb

    re = N.cos(bb.real*l)
    im = -N.sin(bb.real*l)
    phi = phase(re + 1j*im, adjust_range, perform_mod) 

    return phi
    
lagrange_order = 4
mesh_refine_factor = 6
mesh = Rectangle(0, 0, l, a, int(N.ceil(l/a)), 1)
for r in range(mesh_refine_factor):
    mesh = refine(mesh)
    pass

# Function space used for constant functions
DG = FunctionSpace(mesh, "DG", 0)
DG_1 = FunctionSpace(mesh, "DG", 1)
# Function space used for solution
V = FunctionSpace(mesh, "Lagrange", lagrange_order)
v = TestFunction(V)
w = TrialFunction(V)

#Number of modes to use to represent the field at each port
M = 3

#number of ports
num_ports = 2
#Define the ports and mark the subdomains
ports = MeshFunction('uint', mesh, 1)
port1 = InputPort()
port2 = OutputPort()

port1.mark(ports, 1)
port2.mark(ports, 2)

#Define the PEC boundary
pec = WaveguideWalls()
zero = Expression("0.0")
bc = DirichletBC(V, zero, pec)

#Define the sinusoidal part of the excitation function
mode_functions = {}
for m in range(M):
    mode_functions[m] = Expression('sin(m*pi/a*x[1])', {'pi':N.pi, 'm':m+1, 'a':a})



#Define a function for the cutoff squared
#k_o_squared = CutoffSquared(DG)
#k_o_squared.set_frequency(0.0)
k_o_squared = Expression("value", {"value" : 0.0})

c_m = {}


C = N.zeros((V.dim(), num_ports*M))
D = N.zeros((num_ports*M, V.dim()))
coefficients_C = N.zeros(num_ports*M, dtype=N.complex128)

# Assemble C, H, and D - invariant over frequency range
for m in range(M):
    c_m = v*mode_functions[m]
    C_m1 = assemble(c_m*ds(1), exterior_facet_domains=ports)
    C_m2 = assemble(c_m*ds(2), exterior_facet_domains=ports)
    if m == 0:
        H_inc = C_m1.copy()

    d_m = mode_functions[m]*w
    D_m1 = assemble(d_m*ds(1), exterior_facet_domains=ports)
    D_m2 = assemble(d_m*ds(2), exterior_facet_domains=ports)

#    set the entries of the C and D matrices
    k = 0
    C[:,k*M + m] = C_m1.array()[:]
    D[k*M + m, :] = D_m1.array()[:]
    k = 1
    C[:,k*M + m] = C_m2.array()[:]
    D[k*M + m, :] = D_m2.array()[:]


# Define the form for f
s_ii = dot(grad(v), grad(w))*dx
t_ii = k_o_squared*v*w*dx

f_ii = s_ii - t_ii

S = assemble(s_ii)
T = assemble(t_ii)

# creat the LHS and RHS matrices
LHS = scipy.sparse.csr_matrix ((V.dim()+num_ports*M, V.dim()+num_ports*M), dtype=N.complex128 )
RHS = N.zeros((V.dim()+num_ports*M, 1), dtype=N.complex128)

# the entries corresponding to D are invariant under frequency
#LHS[:num_ports*M,num_ports*M:] = D
AD = scipy.sparse.hstack( (scipy.sparse.eye(M*num_ports, M*num_ports, 0, N.complex128, 'csr' ), scipy.sparse.csr_matrix ( D, dtype=N.complex128 ) ), 'csr' )

f_range = [10., 11., 12.,]
#f_range = N.arange(10.0,15.25,0.25).tolist()
#f_range = [1.0,2.0]
#f_range = N.arange(10.0,15.0,0.5).tolist()
f_scale = 1e6

results = N.zeros((1+num_ports*M,len(f_range)), dtype=N.complex128)
F = uBLASSparseMatrix()

for i in range(len(f_range)):
    f = f_range[i]*f_scale
#    k_o_squared.set_frequency(f)
#    k_o_squared.f = f
    k_o_squared.value = k_o(f)**2
    
    print "Assemble F"
#    F = S - k_o_squared*T
    assemble(f_ii, tensor=F )
#    F = assemble(s_ii - k_o(f)**2*t_ii)
    print "apply BC"
    bc.apply(F)
    print "Done"
    
#    P.show()
        
    print "Setting matrix entries"
    t = Timer("Set matrix entries")
    Beta = Beta_km(f)
    ko = k_o(f)
    coefficient_H_C = 1j*Beta*2/N.sqrt(a*b)*N.sqrt(ko*Z_0/Beta)
    coefficient_A_E = N.sqrt(a/b)*N.sqrt(ko*Z_0/Beta)
    
    for k in range(num_ports):
        for m in range(M):
            j = k*M + m
            
            B_km = Beta_km(f, m+1)
            AD[j,j] = -N.sqrt(a/b)*N.sqrt(ko*Z_0/B_km)
            coefficients_C[j] = 1j*B_km*2/N.sqrt(a*b)*N.sqrt(ko*Z_0/B_km)
            
    
    F_sp = dolfin_ublassparse_to_scipy_csr ( F ) # scipy.sparse.csr_matrix ( F.array() )
    C_sp = scipy.sparse.csr_matrix ( N.dot(C, N.diag(coefficients_C)) )

#    LHS[num_ports*M:,num_ports*M:] = F.array()[:,:]
#    LHS[num_ports*M:,:num_ports*M] = N.dot(C, N.diag(coefficients_C))
    LHS = scipy.sparse.vstack( ( AD, scipy.sparse.hstack( (C_sp, F_sp ), 'csr') ), 'csr' )
    
    RHS[0,0] = -LHS[0,0]
    RHS[num_ports*M:,0] = coefficient_H_C*H_inc.array()[:]
    t.stop()
    
    print "Solving (sparse)"
    t = Timer("Sparse Solve")
    BE = solve_sparse_system ( LHS, RHS )
    t.stop()
    print "Done"
    results[0,i] = f
    
    for k in range(num_ports):
        for m in range(M):
            results[1+k*M+m,i] = BE[k*M+m] 
    
#    print BE[:num_ports*M]

S11 = results[1+0*M,:]
S21 = results[1+1*M,:]

for i in range(results.shape[1]):
    print "f = ", results[0,i]
#    print "|R| = %f <R = %f (%f + j%f)" % (abs(results[1,i]), phase(results[1,i]), results[1,i].real, results[1,i].imag)
#    print "|T| = %f <T = %f (%f + j%f)" % (abs(results[2,i]), phase(results[2,i]), results[2,i].real, results[2,i].imag)
#    print "|R| + |T| = %f" % (abs(results[1,i]) + abs(results[2,i]))
    print "|R| = %f <R = %f (%f + j%f)" % (abs(S11[i]), phase(S11[i]), S11[i].real, S11[i].imag)
    print "|T| = %f <T = %f (%f + j%f)" % (abs(S21[i]), phase(S21[i]), S21[i].real, S21[i].imag)
    print "|R| + |T| = %f" % (abs(S11[i]) + abs(S21[i]))

print RHS.shape
summary()


#S11 = N.zeros_like(results[1+0*M,:])
#S21 = N.zeros_like(results[1+1*M,:])
#
#for m in range(M):
#    S11 += results[1+0*M+m,:]
#    S21 += results[1+1*M+m,:]
    



phase_ref = reference_phase(f_range, False, perform_mod=False)
plot_two_axis_data(results[0,:].real/1e6, [N.zeros(results[0,:].shape), magnitude(S21)], 
                   [phase_ref, phase(S21, False, perform_mod=False),],
                   Y_labels = ['Magnitude [dB]', '  Phase [degrees]'], 
                   Y1_styles = ['ks', '-', ], 
                   Y1_limits = [-20, 1],
                   Y2_styles = ['ko', '--',])

P.xlabel('Frequency [GHz]')
#P.savefig('results/images/hollow_h_plane_S_parameters.eps')
P.show()


# Analytical

#
#mag = Function(V)
#mag.vector().set(N.abs(BE[2:,0]))
#
#plot(mag)
#
#tmp = N.arctan(N.imag(BE[2:,0])/N.real(BE[2:,0]))*180/N.pi
#
#phase = Function(V)
#phase.vector().set(tmp)
#
#plot(phase, rescale=True)
#
#interactive()
