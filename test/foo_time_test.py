import time
import numpy as np
import scipy as sp
from scipy import sparse
from scipy.sparse.linalg import eigs, eigsh

def main(aff_arr):
    n = 8100
    v = np.random.rand(n, 1)
    aff_arr_rowsum = aff_arr.sum(axis=0)
    D_sqrt_inv = sp.sparse.spdiags(np.reciprocal(np.sqrt(aff_arr_rowsum)), 0, n, n)


    tic = time.time()
    for i in range(0, 300):
        Av1 = D_sqrt_inv.dot(v)
        Av2 = aff_arr.dot(Av1)
        Av3 = D_sqrt_inv.dot(Av2)
    toc = time.time()
    print(toc - tic)

    tic = time.time()
    M = D_sqrt_inv*aff_arr*D_sqrt_inv
    toc = time.time()
    print(toc - tic)


    tic = time.time()
    [d1, V1] = eigs(M, k = 1, which='LR', tol = 1e-7)
    toc = time.time()
    print(toc - tic)


    tic = time.time()
    [d2, V2] = eigsh(M, k = 1, tol = 1e-7)
    toc = time.time()
    print(toc - tic)


    return d1, V1, d2, V2


