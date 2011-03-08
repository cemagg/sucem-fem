"""
Test driver code. Anything of lasting importance should be factored out pronto
u"""
from __future__ import division
from itertools import izip
import os, sys
import pickle
import random
import numpy as N
import scipy
import gc
import time as t
from scipy.io import mmread
from NewCode.Utilities import Struct

precond_type = sys.argv[1]
filename = sys.argv[2]
print "Loading matrix"
A = mmread(file(filename+'.mtx'))
print "Done"
A_arrays = pickle.load(file(filename+'-arrays.pickle'))
RHS = A_arrays.RHS; imp_len = A_arrays.imp_len
cur_val = A_arrays.dofarrays[0][0:imp_len]
next_val = A_arrays.dofarrays[1][0:imp_len] # For comparison

from scipy.sparse.linalg.dsolve import factorized


precon_opts = Struct(
    cg=Struct(ksp_type='cg',
              pc_type='none'),
    cg_icc0=Struct(ksp_type='cg',
                   pc_type='icc'),
    cg_icc0pd=Struct(ksp_type='cg',
                     pc_type='icc',
                     pc_factor_shift_positive_definite=1),
    cg_icc1=Struct(ksp_type='cg',
                   pc_type='icc',
                   pc_factor_levels=1),
    cg_icc1pd=Struct(ksp_type='cg',
                     pc_type='icc',
                     pc_factor_shift_positive_definite=1,
                     pc_factor_levels=1),
    cg_icc2=Struct(ksp_type='cg',
                   pc_type='icc',
                   pc_factor_levels=2),
    cg_icc2pd=Struct(ksp_type='cg',
                     pc_type='icc',
                     pc_factor_shift_positive_definite=1,
                     pc_factor_levels=2),
    cg_icc3=Struct(ksp_type='cg',
                   pc_type='icc',
                   pc_factor_levels=3),
    cg_icc3pd=Struct(ksp_type='cg',
                     pc_type='icc',
                     pc_factor_levels=3,
                     pc_factor_shift_positive_definite=1),
    cg_amg=Struct(ksp_type='cg',
                  pc_type='hypre',
                  pc_hypre_type='boomeramg'),
    cg_spai=Struct(ksp_type='cg',
                  pc_type='spai',)
                  )

#     PETSC_DECIDE      =   -1
#     OptDB = PETSc.Options()
#     for key in OptDB.getAll().keys(): OptDB.delValue(key)
#     OptDB[''] = 'cg'
#     # OptDB['ksp_type']          = 'gmres'
#     # OptDB['ksp_gmres_restart'] = 50
#     # OptDB['pc_type'] = 'none'
#     #OptDB['pc_type']  = 'icc'
#     #OptDB['pc_type']  = 'ilu'
#     #OptDB['pc_factor_levels'] = 5
#     #OptDB['pc_factor_shift_positive_definite'] = 1
#     #OptDB['pc_factor_shift_nonzero'] = PETSC_DECIDE
#     OptDB['pc_type']        = 'hypre'
#     OptDB['pc_hypre_type']  = 'boomeramg'
#     OptDB['pc_hypre_boomeramg_strong_threshold'] = 0.5
#     OptDB['pc_hypre_boomeramg_max_iter'] = 2

if precond_type == 'umf':
    A_csc = A.tocsc()
    gc.collect()
    t.sleep(1)
    t1 = t.time()
    A_lu = factorized(A_csc)
    t2 = t.time()
    x_umf = A_lu(RHS)
    t3 = t.time()
    x_umf = A_lu(RHS)
    t4 = t.time()
    umfcon = A_lu.func_closure[1].cell_contents
    err = scipy.linalg.norm(A.matvec(N.array(x_umf)) - RHS)
    print 'LU decomp time: %f s' % (t2-t1)
    print 'Solve time 1  : %f s' % (t3-t2)
    print 'Solve time 2  : %f s' % (t4-t3)
    print 'Err           : %f' % err
    print 'dofs: %d  nnz: %d' % (A.shape[0], A.nnz)
else:
    guess = int(sys.argv[3])
    A = A.tocsr()
    A.ensure_sorted_indices(inplace=True)
    import petsc4py
    petsc4py.init(sys.argv)
    from petsc4py import PETSc
    OptDB = PETSc.Options()
    for k,v in precon_opts[precond_type].items():
        OptDB[k]=v
    A_petsc = PETSc.Mat().createAIJ(A.shape,
                                    csr=(A.indptr,
                                         A.indices,
                                         A.data))
    # obtain vectors for storing
    # the solution  and the rhs
    x_petsc, b_petsc = A_petsc.getVecs()
    # fill the rhs PETSc vector
    # from the rhs numpy array
    b_petsc[...] = RHS # note the syntax sugar
    ksp = PETSc.KSP().create()
    ksp.setFromOptions()  # configure from OptDB
    
    ksp.setInitialGuessNonzero(guess)
    ksp.setOperators(A_petsc)   
    ksp.setTolerances(1e-10)
    x_petsc[...] = cur_val
    ksp.solve(b_petsc, x_petsc) 
    x_petsc[...] = cur_val
    gc.collect()
    t.sleep(1)
    t1 = t.time()
    ksp.solve(b_petsc, x_petsc) 
    t2 = t.time()
    x_petsc[...] = cur_val
    t3 = t.time()
    ksp.solve(b_petsc, x_petsc)
    t4 = t.time()
    x_iter = x_petsc[...]
    err = scipy.linalg.norm(A.matvec(N.array(x_iter)) - RHS)
    print 'Solve time    : %f s' % (t2-t1)
    print 'Solve time    : %f s' % (t4-t2)
    print 'Err           : %f' % err
    print ksp.getPC().view()
    print 'dofs: %d  nnz: %d' % (A.shape[0], A.nnz)


