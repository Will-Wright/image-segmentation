# This function verifies the function SolveSDPwithCvxopt
# is correctly transforming SDPs from our form to
# the cvxopt form.

import numpy as np
from scipy.linalg import orth
import sys
import os
sys.path.append(os.getcwd() + "/src")
sys.path.append(os.getcwd() + "/src/SDPSubspaceSolver")
sys.path.append(os.getcwd() + "/src/SDPSubspaceSolver/util")
sys.path.append(os.getcwd() + "/src/NewtonEigSolver")
from SDPSubspaceSolver import *
from SolveSDPwithCvxopt import *
from ComputeAforward import *
from ComputeAadjoint import *
from ProjectModelOntoSubspace import *
from NewtonEigSolver import *
from scipy.sparse.linalg import eigsh
import time


def main(n=300, m=1, cvx_verbose=False, cvx_tol_scale=1, \
         run_cvx_on_full_SDP = True, seed_num=0, \
         eig_solver = "lobpcg", eig_solver_Newton = "eigsh", \
         pr_dim_update=[], run_Newton_solver = False, use_Newton_step = False):
#    A = np.zeros((3, 3, 2))
#    A[:,:, 0] = np.array( [ [2., 6., 4.], [6., -2., 7.], [4., 7., 10.]  ] )
#    A[:,:, 1] = np.array( [ [1., 7., 6.], [7., 14., 9.], [6., 9., 8.]  ] )
#    b = np.array([[9.], [15.]])
#    C = np.array( [ [5., 2., 1.], [2., 7., 3.], [1., 3., 8.] ] )
#    X = np.array( [ [0.152513984269634, 0.206811764229985, 0.175877250073743],
#            [0.206811764229985, 0.280440581743945, 0.238492787099002],
#            [0.175877250073743, 0.238492787099002, 0.202819500417786] ] )
#    y = np.array( [ [0.183267184630601], [0.353917564707084] ] )
#    Z = np.array( [ [4.279548066242836, -1.577026060733195, -1.856574126764909],
#            [-1.577026060733195, 2.411688463573145, -1.468128374777965],
#            [-1.856574126764909, -1.468128374777965, 3.335987636248442] ] )

    np.random.seed(seed_num)

    print("Generating random test SDP")
    C = get_rand_spdmat(n, 0, seed_num, is_PSD=True)
    b = np.sqrt(n)*np.random.rand(m,1)

    if run_Newton_solver:
        b[0] = 1
        b[1] = 0

    A_arr = np.zeros((n, n, m))
    for j in range(0, m):
        A_arr[:,:,j] = get_rand_spdmat(n, j, seed_num+1+j, \
                                       is_PSD=False, \
                                       run_Newton_solver = run_Newton_solver)
    A = (A_arr[:,:,0],)
    for j in range(1, m):
        A = A + (A_arr[:,:,j],)

    if run_cvx_on_full_SDP:
        print('Solving test SDP with cvxopt\n')
        time_start = time.time()
        (X, y) = SolveSDPwithCvxopt(A_arr, b, C, verbose=True, \
                                    gap_tol=1e-7*cvx_tol_scale, \
                                    rel_tol=1e-6*cvx_tol_scale, \
                                    feas_tol=1e-7*cvx_tol_scale)
        runtime_cvx = time.time() - time_start
        print("Rank(X) = ", np.linalg.matrix_rank(X, tol=1e-7))
    else:
        X = []
        y = []

    print('Solving test SDP with SDPSubspaceSolver\n')
    (X_SDPSS_sub, V_SDPSS_sub, y_SDPSS, SDPSS_data) \
    = SDPSubspaceSolver(A, b, C, eig_solver = eig_solver, \
                        eigsh_tol_init=1e-8, eigsh_maxiters = 20000, \
                        eigsh_tol_update=1e-2, eigsh_tol_min=1e-14, \
                        lobpcg_tol_init=1e-8, lobpcg_tol_update=1e-2, \
                        lobpcg_tol_min=1e-14, lobpcg_maxiters = 100000, \
                        pr_dim_update=pr_dim_update, \
                        cvx_verbose=cvx_verbose, return_data=True)

    if run_Newton_solver and (m==2):
        (V1, f, NES_data) \
        = NewtonEigSolver(A, b, C, eig_solver = eig_solver_Newton, \
                       eigsh_tol_init=1e-8, eigsh_maxiters = 20000, \
                       eigsh_tol_update=1e-2, eigsh_tol_min=1e-14, \
                       lobpcg_tol_init=1e-8, lobpcg_tol_update=1e-2, \
                       lobpcg_tol_min=1e-14, lobpcg_maxiters = 100000, \
                       use_Newton_step=use_Newton_step, return_data=True)

    print('Computing residuals for experiments\n')

    # Computes residuals for SDPSS
    (A_sub, C_sub) = ProjectModelOntoSubspace(A, C, V_SDPSS_sub)
    pr_res_SDPSS = np.linalg.norm(ComputeAforward(A_sub,X_SDPSS_sub) - b)
    Aty_SDPSS = ComputeAadjoint(A,y_SDPSS)
    [du_res_SDPSS, V] = eigsh(C-Aty_SDPSS, k = m, which = 'SA', tol = 1e-14, maxiter=20000)
    pr_obj = np.tensordot(C_sub,X_SDPSS_sub,2)
    du_obj = np.tensordot(b,y_SDPSS,2)
    gap_res_SDPSS = pr_obj - du_obj # must be >= 0 if pr and du feasible

    print('         | primal residual |  dual residual  |  duality gap   | runtime secs')
    print('  SDPSS ', '|    %1.3e   '%pr_res_SDPSS, '|   % 1.3e   '%du_res_SDPSS[0], \
          '|   %1.3e' %gap_res_SDPSS, '   |  %11.2f' %SDPSS_data['time_total'])
    if run_Newton_solver and (m==2):
        print(' Newton ', '|    %1.3e   '%pr_res_SDPSS, '|   % 1.3e   '%du_res_SDPSS[0], '|   %1.3e' %gap_res_SDPSS, '   |  %11.2f' %SDPSS_data['time_total'])
    if run_cvx_on_full_SDP:
        pr_res = np.linalg.norm(ComputeAforward(A,X) - b)
        Aty = ComputeAadjoint(A,y)
        [du_res, V] = eigsh(C-Aty, k = 3, which = 'SA', tol = 1e-14, maxiter=20000)
        pr_obj = np.tensordot(C,X,2)
        du_obj = np.tensordot(b,y,2)
        gap_res = pr_obj - du_obj # must be >= 0 if pr and du feasible
        print('full SDP', '|    %1.3e   '%pr_res, '|   % 1.3e   '%min(du_res), \
              '|   %1.3e' %gap_res, '   |  %11.2f' %runtime_cvx)

    return X, y, A, b, C, X_SDPSS_sub, V_SDPSS_sub, y_SDPSS


def get_rand_spdmat(n, j, seed_num, is_PSD=True, run_Newton_solver = False):
    np.random.seed(seed_num)

    U = orth(np.random.rand(n,n))
    if run_Newton_solver and j == 0:
        D = np.ones((n)) + np.diag(np.random.rand(n))
    elif run_Newton_solver and j == 1:
        D = -0.5*np.ones((n)) + np.diag(np.random.rand(n))
    else:
        D = np.ones((n)) + np.diag(np.random.rand(n))

    M = np.dot(np.dot(U,D), U.transpose())
    M = M + M.transpose()
    return M
