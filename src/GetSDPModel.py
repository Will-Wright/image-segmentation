# Computes semidefinite programming model for image segmentation problem

import numpy as np
import scipy as sp
from scipy.sparse import spdiags
from scipy.linalg import orth
from scipy.sparse.linalg import eigsh

def GetSDPModel(aff_arr, N, eig_tol=1e-4):

    use_NCuts_vecs = False
    num_NCuts_vecs = 1
    use_minmax_vecs = True

    n = aff_arr.shape[0]
    aff_arr_rowsum = aff_arr.sum(axis=0)
    D = sp.sparse.spdiags(aff_arr_rowsum, 0, n, n)
    DmW = D - aff_arr

    C = N.transpose()*DmW*N
    b = np.array([[1.0],[0.0]])
    A1 = N.transpose()*D*N

    a2_arr = np.concatenate( ((-1/n)*np.ones(n-1,), np.array([(n-1)/n]) ) )
    A2_full = sp.sparse.spdiags(a2_arr, 0, n, n)
    A2 = N.transpose()*A2_full*N

    A = (A1, A2)

    N_cols = N.shape[1]
    V_init = N[:,N_cols-1].toarray()

    if use_NCuts_vecs:
        aff_arr_rowsum_sub = aff_arr_rowsum[0:n-1]
        aff_arr_sub = aff_arr[0:n-1, 0:n-1]
        V_minmax = np.zeros((n, num_NCuts_vecs))
        D_sub_sqrt_inv \
          = sp.sparse.spdiags(np.reciprocal(np.sqrt(aff_arr_rowsum_sub)), 0, n-1, n-1)
        A_NCuts_sub = D_sub_sqrt_inv*aff_arr_sub*D_sub_sqrt_inv
        [d, V_minmax[0:n-1,0:num_NCuts_vecs]] = eigsh(A_NCuts_sub,\
                                                      k = num_NCuts_vecs,\
                                                      which = 'LA', tol = eig_tol)
        V_init = orth(np.concatenate((V_init, V_minmax), axis=1))

    if use_minmax_vecs:
        if 0==aff_arr_rowsum[0,n-1]:
            aff_arr_rowsum[0,n-1] = 1e16
        D_inv = sp.sparse.spdiags(np.reciprocal(aff_arr_rowsum), 0, n, n)
        [d, v_min_sub] = eigsh(D_inv, k = 1, which = 'SA', tol = eig_tol)
        [d, v_max_sub] = eigsh(D_inv, k = 1, which = 'LA', tol = eig_tol)
        v_min = N.transpose()*v_min_sub
        v_max = N.transpose()*v_max_sub
        V_init = orth(np.concatenate((V_init, v_min_sub, v_max_sub), axis=1))

    V_N_init = N.transpose()*V_init

    return A, b, C, V_N_init
