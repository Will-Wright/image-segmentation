# This function solves the following eigenvalue problem
# using a Newton method.
#
#         x'(C + t*A2)x
# max min -------------
#  t   x      x'A1x
#
# A is a linear operator defined by
# 	forward: A(X) = vec( tr(A1'*x), ..., tr(Am'*X) )
# 	adjoint: A'(y) = y1*A1 + ... + ym*Am
#	Ai are all n-by-n symmetric or Hermitian
# C is a symmetric or Hermitian matrix
# b is a real m-vector
# X is the primal variable, an n-by-n symmetric or Hermitian matrix
# y is the dual variable, a real m-vector

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

from GetPrimalSubspace import *


def NewtonEigSolver(A, b, C, V_init=np.zeros((0,0)), max_iters=50, \
                    verbosity=1, \
                    pr_tol = 1e-8, du_tol = 1e-8, gap_tol = 1e-6, \
                    eig_solver = "lobpcg", use_Newton_step = False, \
                    eigsh_tol_init=1e-4, eigsh_tol_update=1e-2, eigsh_tol_min=1e-10, \
                    eigsh_maxiters = 100000, lobpcg_tol_init=1e-5, \
                    lobpcg_tol_update=1e-2, lobpcg_tol_min=1e-14, \
                    lobpcg_maxiters = 100000, return_data=False):
    time_start = time.time()

    if sp.sparse.issparse(A[0]) or sp.sparse.issparse(A[1]) \
        or sp.sparse.issparse(C):
        is_sparse = True
    else:
        is_sparse = False

    time_eig = 0
    time_lin_solve = 0
    time_main = 0
    time_total = 0
    num_eig_calls = 0
    num_lin_solves = 0
    n = C.shape[0]
    m = len(b)
    # Initialize parameters
    eigsh_tol = eigsh_tol_init
    lobpcg_tol= lobpcg_tol_init

    SDP_pr_tol = pr_tol
    SDP_du_tol = du_tol
    SDP_gap_tol = gap_tol

    converged = False

    attempt_Newton_step = False
    take_Newton_step = False

    if return_data:
        V_data = list()
        y_data = list()
    if V_init.shape[1] > 0:
        V_data.append(V_init)


    # Initial iteration
    iter = 0


    time_main = time_main + time.time() - time_start; time_start = time.time()
    # Initialize subspace (attempt to find subspace where A(X_sub) = b is feasible)
    v0 = np.random.rand(n,1)
    v0 = v0 / np.linalg.norm(v0)
    (a, V1) = GetPrimalSubspace(A=A[1], B=A[0], num_eigs=1, which = 'SA', \
                                eigsh_tol=eigsh_tol, eigsh_maxiters=eigsh_maxiters,
                                lobpcg_tol=lobpcg_tol,
                                lobpcg_maxiters=lobpcg_maxiters, \
                                eig_solver = "eigsh", v_prev = v0)
    num_eig_calls = num_eig_calls + 1
    if return_data:
        V_data.append(V1)
    v0 = np.random.rand(n,1)
    v0 = v0 / np.linalg.norm(v0)
    (c, Vn) = GetPrimalSubspace(A=A[1], B=A[0], num_eigs=1, which = 'LA', \
                                eigsh_tol=eigsh_tol, eigsh_maxiters=eigsh_maxiters,
                                lobpcg_tol=lobpcg_tol, \
                                lobpcg_maxiters=lobpcg_maxiters, \
                                eig_solver = "eigsh", v_prev = v0)
    num_eig_calls = num_eig_calls + 1
    if return_data:
        V_data.append(Vn)
    time_eig = time_eig + time.time() - time_start; time_start = time.time()

    vCv = V1.transpose().dot(C.dot(V1))
    vA1v = V1.transpose().dot(A[0].dot(V1))
    vA2v = V1.transpose().dot(A[1].dot(V1))
    vnA1vn = Vn.transpose().dot(A[0].dot(Vn))

    b = vCv / (vA1v + 1e-8)
    d = Vn.transpose().dot(C.dot(Vn)) / (vnA1vn + 1e-8)

    if verbosity==2:
        print(a,b,c,d)

    t = (d-b)/(a-c);
    f = c*t+d;
    V1 = V1 + Vn
    V1 = V1 / np.linalg.norm(V1)
    vCv = V1.transpose().dot(C.dot(V1))
    vA1v = V1.transpose().dot(A[0].dot(V1))
    vA2v = V1.transpose().dot(A[1].dot(V1))

    # prevents Newton step until 1st order step
    # taken and proper values found for t_min and t_max
    t_min = 1e16; t_max = -1e16;
    t_diff = abs(t_min - t_max);
    g_min = a; k_min = b;
    g_max = c; k_max = d;

    time_main = time_main + time.time() - time_start; time_start = time.time()
    if verbosity>=1:
        if verbosity==2:
            print('NewtonSolver called with eigenvalue solver', eig_solver)
            print('    |                                N-step | Runtime Ratios |                 | ')
            print(' it |  probj       prres    dugap   try take| eig  Ax=b main | eigtol     time |  f            t')
        else:
            print('NewtonSolver called with eigenvalue solver', eig_solver)
            print('    |                                N-step | Runtime Ratios | ')
            print(' it |  probj       prres    dugap   try take| eig  Ax=b main | eigtol     time')

    time_main = time_main + time.time() - time_start;
    # Skips runtime for following eigenvalue problem,
    # not in original method, just used for convergence msmt
    du_gap = np.absolute(vCv[0] - f)
    pr_res = np.linalg.norm([vA1v-1, vA2v-0])

    if is_sparse:
        fA1 = A[0].multiply(f)
        tA2 = A[1].multiply(t)
    else:
        fA1 = f*A[0]
        tA2 = t*A[1]

#    (du_res, V_du_res) \
#        = GetPrimalSubspace(A=C - fA1 + tA2, \
#                            num_eigs=1, which = 'SA', \
#                            eigsh_tol=1e-10, eigsh_maxiters=100000, \
#                            lobpcg_tol=1e-10, lobpcg_maxiters=100000, \
#                            eig_solver = "eigsh", v_prev = V1)
    if (pr_res<=SDP_pr_tol) and (du_gap[0]<=SDP_gap_tol):
        converged = True
    time_start = time.time()

    time_iter = time_eig+time_main+time_lin_solve
    time_total = time_total + time_iter

    if eig_solver == "eigsh":
        tol_print = eigsh_tol
    else:
        tol_print = lobpcg_tol
    if verbosity>=1:
        if verbosity==2:
            print("%3i" % iter, \
                  "| % 3.5e" %f, \
        #          "% 3.2e"%du_res[0], \
                  "%3.2e"%pr_res, \
                  "%3.2e"%du_gap[0], \
                  "% 2s " % "N", \
                  "% 2s |" % "N", \
                  "%4.2f"%(time_eig/ (time_iter) ), \
                  "%4.2f"%(time_lin_solve/ (time_iter) ), \
                  "%4.2f"%(time_main/ (time_iter) ), \
                  "| %3.1e" %tol_print, \
                  " %6.2f"% (time_iter), \
                  "| % 3.5e"% f, \
                  "% 3.5e" % t )
        else:
            print("%3i" % iter, \
                  "| % 3.5e" %f, \
        #          "% 3.2e"%du_res[0], \
                  "%3.2e"%pr_res, \
                  "%3.2e"%du_gap[0], \
                  "% 2s " % "N", \
                  "% 2s |" % "N", \
                  "%4.2f"%(time_eig/ (time_iter) ), \
                  "%4.2f"%(time_lin_solve/ (time_iter) ), \
                  "%4.2f"%(time_main/ (time_iter) ), \
                  "| %3.1e" %tol_print, \
                  " %6.2f"% (time_iter))

    if return_data:
        y_data.append( (f,t) )

    iter = iter + 1
    time_main = 0; time_eig = 0; time_lin_solve = 0;

    # Main loop
    while (iter <= max_iters) and (converged == False):

        eigsh_tol = max(eigsh_tol_update*eigsh_tol, eigsh_tol_min)
        lobpcg_tol = max(lobpcg_tol_update*lobpcg_tol, lobpcg_tol_min)
        time_main = time_main + time.time() - time_start; time_start = time.time()

        if is_sparse:
            tA2 = A[1].multiply(t)
        else:
            tA2 = t*A[1]
        (f, V1) = GetPrimalSubspace(A=C + tA2, B=A[0], num_eigs=1, which = 'SA', \
                                    eigsh_tol=eigsh_tol, \
                                    eigsh_maxiters=eigsh_maxiters,
                                    lobpcg_tol=lobpcg_tol, \
                                    lobpcg_maxiters=lobpcg_maxiters, \
                                    eig_solver = eig_solver, v_prev = V1)
        num_eig_calls = num_eig_calls + 1
        time_eig = time_eig + time.time() - time_start; time_start = time.time()
        if return_data:
            V_data.append(V1)

        vA1v = V1.transpose().dot(A[0].dot(V1))
        vA2v = V1.transpose().dot(A[1].dot(V1))
        vCv = V1.transpose().dot(C.dot(V1))
        k_new = vCv / vA1v;
        g_new = vA2v / vA1v;

        if g_new < 0:
            k_min = k_new;
            g_min = g_new;
            t_max = t;
        else:
            k_max = k_new;
            g_max = g_new;
            t_min = t;

        t_grad_step = (k_max - k_min) / (g_min - g_max);
        t_prev = t;
        time_main = time_main + time.time() - time_start;

        # Skips runtime for following eigenvalue problem,
        # not in original method, just used for convergence msmt
        du_gap = np.absolute(vCv[0] - f[0])
        pr_res = np.linalg.norm([vA1v-1, vA2v-0])

        if is_sparse:
            fA1 = A[0].multiply(f)
            tA2 = A[1].multiply(t)
        else:
            fA1 = f*A[0]
            tA2 = t*A[1]
#        (du_res, V_du_res) \
#            = GetPrimalSubspace(A=C - fA1 + tA2, \
#                                num_eigs=1, which = 'SA', \
#                                eigsh_tol=1e-10, eigsh_maxiters=100000, \
#                                lobpcg_tol=1e-10, lobpcg_maxiters=100000, \
#                                eig_solver = "eigsh", v_prev = V1)

        if (pr_res<=SDP_pr_tol) and (du_gap[0]<=SDP_gap_tol):
            converged = True
        time_start = time.time()

        if not converged:
            if use_Newton_step and (t_min < 1e16) and (t_max > -1e16):
                attempt_Newton_step = True

                if is_sparse:
                    gA1 = A[0].multiply(g_new)
                    Vt = sp.sparse.csr_matrix(V1.transpose())


                    A_Newton = sp.sparse.vstack( [C - fA1 + tA2, \
                                                  Vt] )
                    b_Newton = sp.sparse.vstack( [-A[1].dot(V1) + (gA1).dot(V1),\
                                                  np.zeros((1,1))] )

                    AtA = A_Newton.transpose().dot(A_Newton)
                    Atb = A_Newton.transpose().dot(b_Newton)
                    v_deriv = sp.sparse.linalg.spsolve(AtA, Atb)

                else:
                    gA1 = (vA2v/vA1v)*A[0]
                    Vt = V1.transpose()
                    A_Newton = np.concatenate( (C - fA1 + tA2, \
                                                Vt \
                                                ), \
                                                axis=0\
                                         )
                    b_Newton = np.concatenate( (-A[1].dot(V1) + (gA1).dot(V1), \
                                                np.zeros((1,1)) \
                                               ),\
                                               axis=0\
                                         )
                    AtA = A_Newton.transpose().dot(A_Newton)
                    Atb = A_Newton.transpose().dot(b_Newton)
                    v_deriv = np.linalg.solve(AtA, Atb)

                W = (A[1] - gA1).dot(v_deriv)
                h = (2/vA1v) * ( V1.transpose().dot( W ) );
                t_Newton = -(g_new / h) + t;
                if (t_Newton >= t_min) and (t_Newton <= t_max):
                    take_Newton_step = True
                    t = t_Newton;
                else:
                    take_Newton_step = False
                    t = t_grad_step;

                num_lin_solves = num_lin_solves + 1
                time_lin_solve = time_lin_solve + time.time() \
                                 - time_start; time_start = time.time()

            else:
                attempt_Newton_step = False
                take_Newton_step = False
                t_Newton = 0.0
                t = t_grad_step;

            t_diff = abs(t-t_prev);
            if take_Newton_step:
                tNs = "Y"
            else:
                tNs = "N"
            if attempt_Newton_step:
                aNs = "Y"
            else:
                aNs = "N"

            if return_data:
                y_data.append( (f,t) )
        else:
            tNs = "N"
            aNs = "N"

        time_iter = time_eig+time_main+time_lin_solve
        time_total = time_total + time_iter

        if verbosity>=1:
            if eig_solver == "eigsh":
                tol_print = eigsh_tol
            else:
                tol_print = lobpcg_tol
            if verbosity==2:
                print("%3i" % iter, \
                      "| % 3.5e" %f, \
        #              "% 3.2e"%du_res[0], \
                      "%3.2e"%pr_res, \
                      "%3.2e"%du_gap[0], \
                      "% 2s " % aNs, \
                      "% 2s |" % tNs, \
                      "%4.2f"%(time_eig/ (time_iter) ), \
                      "%4.2f"%(time_lin_solve/ (time_iter) ), \
                      "%4.2f"%(time_main/ (time_iter) ), \
                      "| %3.1e" %tol_print, \
                      " %6.2f"% (time_iter), \
                      "| % 3.5e"% f, \
                      "% 3.5e" % t)
            else:
                print("%3i" % iter, \
                      "| % 3.5e" %f, \
        #              "% 3.2e"%du_res[0], \
                      "%3.2e"%pr_res, \
                      "%3.2e"%du_gap[0], \
                      "% 2s " % aNs, \
                      "% 2s |" % tNs, \
                      "%4.2f"%(time_eig/ (time_iter) ), \
                      "%4.2f"%(time_lin_solve/ (time_iter) ), \
                      "%4.2f"%(time_main/ (time_iter) ), \
                      "| %3.1e" %tol_print, \
                      " %6.2f"% (time_iter) )

        iter = iter + 1
        time_main = 0; time_eig = 0; time_lin_solve = 0;

    if verbosity>=1:
        if converged:
            print("\nEXIT -- Optimal solution found.\n")
        else:
            print("\nEXIT -- Maximum number of iterations.\n")
        print("Total runtime     :  ", "%1.2f" %time_total)
        print("Num eig calls     :  ", "%1i" %num_eig_calls)
        print("Num Ax=b solves   :  ", "%1i" %num_lin_solves)
        print("\n")

    if return_data:
        data_out = {'time_total':time_total,
                    'num_eig_calls' : num_eig_calls,
                    'V_data':V_data,
                    'y_data':y_data}
        out = (V1, f, data_out)
    else:
        out = (V1, f)

    return out


