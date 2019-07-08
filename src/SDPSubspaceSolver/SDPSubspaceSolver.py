# This function solves semidefinite programs (SDPs)
# using a subspace method on the primal variable X.
#
# This function uses the following SDP form with primal
#           min   tr(C'*X)
# (SDP-P)    X
#            st   A(X) = b, X >= 0
#
# and dual
#           max   b'*y
# (SDP-D)    y
#            st   Z = C - A'(y) >= 0
#
# A is a linear operator defined by
# 	forward: A(X) = vec( tr(A1'*x), ..., tr(Am'*X) )
# 	adjoint: A'(y) = y1*A1 + ... + ym*Am
#	Ai are all n-by-n symmetric or Hermitian
# C is a symmetric or Hermitian matrix
# b is a real m-vector
# X is the primal variable, an n-by-n symmetric or Hermitian matrix
# y is the dual variable, a real m-vector
#
# Note: this subspace method is particularly efficient for
# SDPs with a very large objective matrix (e.g., n > 1000)
# and very few primal constraints (e.g., m < 50).
# Well-suited problems include trace-ratio problems, image segmentation
# with normalized cuts, and trust-region subproblems.
#
#
# The following papers are notable in developing this subspace
# method and its theory.
#
# 2002 Oliveira Stewart Soma, A Subspace Semidefinite Programming
# for Spectral Graph Partitioning
# 	- discusses this subspace method for first time in literature
#	- presents a few small applications, no theory
#
# 2015 Kangal Meerbergen Mengi Michiels, A Subspace Method
# for Large Scale Eigenvalue Optimization
#	- reintroduces subspace SDP method as eigenvalue method
#	- establishes basic convergence theory for general case
#
# 2017 Kressner Lu Vandereycken, Subspace acceleration
# for the Crawford number and related eigenvalue optimization problems
#	- establishes super-quadratic (1 + sqrt(2)) local convergence
#  	  for subspace method with m = 1 (only 1 primal constraint)
#
# The following convergence criteria are used for the subspace method
#
#    KKT Condition              Algorithm Criteria
#
#          A(X)  = b     ||A(X) - b||               <= params.primal_residual_tol
#             X >= 0     None: condition holds for all iterates X = Vsub'*Xsub*Vsub
#       C-At(y) >= 0     eig_min(C-At(y))           >= -params.dual_residual_tol
#  <X, C-At(y)>  = 0     None: condition holds for all iterates
#  <C, X> - b'*y = 0     abs(primal_obj - dual_obj) <= params.gap_tol

import numpy as np
import scipy as sp
from scipy.sparse import spdiags
import os, sys, inspect
import time

# Adds cwd to path
cwd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cwd_folder not in sys.path:
    sys.path.insert(0, cwd_folder)
# Adds util subdirectory to path
cwd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.    currentframe() ))[0], 'util')))
if cwd_subfolder not in sys.path:
    sys.path.insert(0, cwd_subfolder)

from PrintBanner import *
from GetPrimalSubspace import *
from ProjectModelOntoSubspace import *
from SolveSDPwithCvxopt import *
from ComputeAadjoint import *
from ComputeAforward import *
from ProjectOntoPrimalSpace import *
from DetermineIfConverged import *
from PrintIteration import *
from UpdateSubspace import *

def SDPSubspaceSolver(A, b, C, V_init=np.zeros((0,0)), \
                      verbosity=1, \
                      max_iters=20, eig_solver = "lobpcg", \
                      pr_tol=1e-8, du_tol=1e-8, gap_tol=1e-6, \
                      eigsh_tol_init=1e-4, eigsh_tol_update=1e-2, \
                      eigsh_tol_min=1e-10, \
                      eigsh_maxiters = 100000, lobpcg_tol_init=1e-5, \
                      lobpcg_tol_update=1e-2, lobpcg_tol_min=1e-14, \
                      lobpcg_maxiters = 100000, \
                      cvx_verbose=False, pr_dim_update=[], \
                      return_data=False, \
                      test_model = False):
    time_start = time.time()

    time_eig = 0
    time_SDP = 0
    time_main = 0
    time_total = 0
    num_eig_calls = 0
    n = C.shape[0]
    m = len(b)
    # Initialize parameters
    eigsh_tol = eigsh_tol_init
    lobpcg_tol= lobpcg_tol_init
    
    #nonsym_err_tol = 1e-8 # tol for warning flag on nonsymmetry of Aty
    converged = False
    pr_dim = V_init.shape[1]
    pr_dim_init = 1
    if (not pr_dim_update):
        pr_dim_update = m
    du_dim = m
    #X = np.zeros((n,n))
    X_sub = np.zeros((pr_dim, pr_dim))
    y = np.zeros((du_dim,1))
    pr_obj = 0
    du_obj = 0
    
    if return_data:
        V_data = list()
        y_data = list()
    if pr_dim > 0:
        V_data.append(V_init)
    
    PrintBanner(eig_solver, verbosity)
    
    # Initial iteration
    iter = 0        
    
    time_main = time_main + time.time() - time_start; time_start = time.time()
                    
    # Initialize subspace (attempt to find subspace where A(X_sub) = b is feasible)        
    if test_model:
        (a, V1) = GetPrimalSubspace(A=A[1], B=A[0], num_eigs=1, which = 'SA', \
                                    eigsh_tol=eigsh_tol, \
                                    eigsh_maxiters=eigsh_maxiters, 
                                    lobpcg_tol=lobpcg_tol, \
                                    lobpcg_maxiters=lobpcg_maxiters, \
                                    eig_solver = eig_solver)
        num_eig_calls = num_eig_calls + 1
        (c, Vn) = GetPrimalSubspace(A=A[1], B=A[0], num_eigs=1, which = 'LA', \
                                    eigsh_tol=eigsh_tol, \
                                    eigsh_maxiters=eigsh_maxiters, 
                                    lobpcg_tol=lobpcg_tol, \
                                    lobpcg_maxiters=lobpcg_maxiters, \
                                    eig_solver = eig_solver)
        num_eig_calls = num_eig_calls + 1
        
        if return_data:
            V_data.append(V1)
            V_data.append(Vn)
        
        if V_sub.shape[1] == 0:
            V_sub = np.concatenate((V1, Vn), axis=1)
        else:
            V_sub = np.concatenate((V_sub, V1, Vn), axis=1)
    elif pr_dim == 0:
        V_sub = V_init
        for j in range(0, m):
            v0 = np.random.rand(n,1)
            v0 = v0 / np.linalg.norm(v0)
            (d, V1_new) = GetPrimalSubspace(A[j], num_eigs=pr_dim_init, \
                                            which = 'SA', \
                                           eigsh_tol=eigsh_tol, \
                                           eigsh_maxiters=eigsh_maxiters, \
                                           lobpcg_tol=lobpcg_tol, \
                                           lobpcg_maxiters=lobpcg_maxiters, \
                                           v_prev = v0, eig_solver = 'eigsh')
            num_eig_calls = num_eig_calls + 1
            if V_sub.shape[1] == 0:
                V_sub = V1_new
            else:
                V_sub = np.concatenate((V_sub, V1_new), axis=1)
            (d, Vn_new) = GetPrimalSubspace(A[j], num_eigs=pr_dim_init, \
                                            which = 'LA', \
                                           eigsh_tol=eigsh_tol, \
                                           eigsh_maxiters=eigsh_maxiters, \
                                           lobpcg_tol=lobpcg_tol,\
                                           lobpcg_maxiters=lobpcg_maxiters, \
                                           v_prev = v0, eig_solver = 'eigsh')
            num_eig_calls = num_eig_calls + 1
            
            if return_data:
                V_data.append(V1_new)
                V_data.append(Vn_new)
                
            V_sub = np.concatenate((V_sub, Vn_new), axis=1)  
    else:
        V_sub = V_init
    time_eig = time_eig + time.time() - time_start; time_start = time.time()

    V_sub = UpdateSubspace(V_sub=V_sub)
    du_res = -1
        
    (A_sub, C_sub) = ProjectModelOntoSubspace(A, C, V_sub)    
    time_main = time_main + time.time() - time_start; time_start = time.time()

    (X_sub, y) = SolveSDPwithCvxopt(A_sub, b, C_sub, verbose=cvx_verbose, gap_tol = 1e-10)
    time_SDP = time_SDP + time.time() - time_start; time_start = time.time()
    
    if return_data:
        y_data.append( y )
    
    Aty = ComputeAadjoint(A,y)
    AX = ComputeAforward(A_sub,X_sub)
    
    pr_obj = np.tensordot(C_sub,X_sub,2)
    pr_res_vec = AX-b
    pr_res = np.linalg.norm(pr_res_vec)
    du_obj = np.tensordot(b,y,2)
    pr_dim = V_sub.shape[1]
    v_prev = np.reshape(V_sub[:,0], (n,1))
    
    converged = DetermineIfConverged(pr_obj, du_obj, pr_res, du_res, pr_tol, du_tol, gap_tol)        
    time_main = time_main + time.time() - time_start; time_start = time.time()
    if eig_solver == "eigsh":
        tol_print = eigsh_tol
    else:
        tol_print = lobpcg_tol
    PrintIteration(iter, pr_obj, du_obj, pr_res, du_res, \
                   time_eig, time_SDP, time_main, pr_dim, tol_print, \
                   y, verbosity)
    time_total = time_total + time_main + time_eig + time_SDP; time_main = 0; time_eig = 0; time_SDP = 0;
    iter = iter + 1
        
#    M = ()

    # Main loop
    while (iter <= max_iters) and (converged == False):
        
        eigsh_tol = max(eigsh_tol_update*eigsh_tol, eigsh_tol_min)
        lobpcg_tol = max(lobpcg_tol_update*lobpcg_tol, lobpcg_tol_min)
        time_main = time_main + time.time() - time_start; time_start = time.time()
#        M = M + (C-Aty,) 
        
        (d, V_new) = GetPrimalSubspace(A=C-Aty, num_eigs=pr_dim_update, \
                                       sigma=None, eigsh_tol=eigsh_tol, \
                                       eigsh_maxiters=eigsh_maxiters, \
                                       lobpcg_tol=lobpcg_tol, \
                                       lobpcg_maxiters=lobpcg_maxiters, \
                                       v_prev = v_prev, eig_solver = eig_solver)
        num_eig_calls = num_eig_calls + 1
        time_eig = time_eig + time.time() - time_start; time_start = time.time()
        du_res = min(d)
        
        if return_data:
            V_data.append(V_new)
        
        #V_new = V_new[:,d<0] # testing indicates best to use all of V_new
        
        V_sub = UpdateSubspace(V_sub=V_sub, V_new=V_new)
        (A_sub, C_sub) = ProjectModelOntoSubspace(A, C, V_sub)
        
        time_main = time_main + time.time() - time_start; time_start = time.time()
        (X_sub, y) = SolveSDPwithCvxopt(A_sub, b, C_sub, verbose=cvx_verbose, gap_tol = 1e-10)
        time_SDP = time_SDP + time.time() - time_start; time_start = time.time()
        
        if return_data:
            y_data.append( y )
        
        Aty = ComputeAadjoint(A,y)
        AX = ComputeAforward(A_sub, X_sub)        
        
        pr_obj = np.tensordot(C_sub, X_sub, 2)                
        pr_res_vec = AX-b        
        pr_res = np.linalg.norm(pr_res_vec)
        du_obj = np.tensordot(b,y,2)
        pr_dim = V_sub.shape[1]
        v_prev = np.reshape(V_new[:,0], (n,1))
        
        converged = DetermineIfConverged(pr_obj, du_obj, pr_res, du_res, pr_tol, du_tol, gap_tol)  
        time_main = time_main + time.time() - time_start; time_start = time.time()
        if eig_solver == "eigsh":
            tol_print = eigsh_tol
        else:
            tol_print = lobpcg_tol        
        PrintIteration(iter, pr_obj, du_obj, pr_res, du_res, \
                       time_eig, time_SDP, time_main, pr_dim, \
                       tol_print, y, verbosity)
        time_total = time_total + time_main + time_eig + time_SDP; time_main = 0; time_eig = 0; time_SDP = 0;
        iter = iter + 1
        
    if verbosity>=1:
        if converged:
            print("\nEXIT -- Optimal solution found.\n")
        else:
            print("\nEXIT -- Maximum number of iterations.\n")
        print("Total runtime   :  ", "%1.2f" %time_total)
        print("Num eig calls   :  ", "%1i" %num_eig_calls)
        print("\n")
        
        
    if return_data:
        data = {'time_total' : time_total,
                'num_eig_calls' : num_eig_calls,
                'V_new' : V_new,
                'V_data' : V_data,
                'y_data' : y_data}
        data_out = (V_data, y_data)
        out = (X_sub, V_sub, y, data)
    else:
        out = (X_sub, V_sub, y)
        
    return out


