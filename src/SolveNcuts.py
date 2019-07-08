# This function solves the normalized cuts problem.
# The function takes as input an affinity matrix W and
# returns the vector corresponding to the smallest
# generalized eigenvalue of the following problem
#
#      x'(D-W)x
# inf -----------   =  eig_min( D-W, D )
#  x     x'Dx
#
# D = Diag(row-sums of W)

import numpy as np
import scipy as sp
from scipy.sparse import spdiags
from scipy.sparse.linalg import eigsh

def SolveNCuts(aff_arr, eig_tol=1e-4):

    invert_D = True
    
    if invert_D:
        # Solves unconstrainted NCuts by inverting D
        n = aff_arr.shape[0]
        aff_arr_rowsum = aff_arr.sum(axis=0)
        D_sqrt_inv = sp.sparse.spdiags(np.reciprocal(np.sqrt(aff_arr_rowsum)), 0, n, n)
        A = D_sqrt_inv*aff_arr*D_sqrt_inv
        [d, V] = eigsh(A, k = 1, which = 'LA', tol = eig_tol)

    else:
    # Solves unconstrained NCuts as generalized eigenvalue problem lam_min(A, B)
        n = aff_arr.shape[0]
        aff_arr_rowsum = aff_arr.sum(axis=0)
        D = sp.sparse.spdiags(aff_arr_rowsum, 0, n, n)
        [d, V] = eigsh(D - aff_arr, M=D, k = 1, which = 'SA', tol = eig_tol)

    print("NCuts routine completed")

    return d, V

