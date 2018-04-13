# This function solves the normalized cuts problem.
# The function takes as input an affinity matrix W and
# returns the vector corresponding to the smallest
# generalized eigenvalue of the following problem
#
#      x'(D-W)x
# inf -----------   =  lam_min( D-W, D )
#  x     x'Dx
#
# D = Diag(row-sums of W)

import numpy as np
import scipy as sp
from scipy.sparse import spdiags
from scipy.sparse.linalg import eigsh

def main(W):

    eig_tol = 1e-4
    [n, m] = W.shape
    d = W.sum(axis=0)
    D = spdiags(d, 0, n, n)

    D_inv = spdiags(1/d, 0, n, n)
    A = D_inv*W

# GOTTA VERIFY matrix A is correct.  not getting A = A' in testing

    A = 0.5*(A + np.transpose(A))

#    n = 1000

#    A = np.random.rand(1000, 1000)
#    A = A + np.transpose(A)
    [d, V] = eigsh(A, k=1, which='SA', tol=eig_tol)
#    D = np.matmul(np.transpose(V),np.matmul(A,V))

    return d, V

