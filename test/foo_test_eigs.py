# Imports necessary packages
import numpy as np
import scipy as sp
import scipy.sparse.linalg

def main():
    A = np.random.rand(8, 8)
    A = np.matmul(np.transpose(A),A)
    [d, V] = sp.sparse.linalg.eigs(A, k = 6)
    D = np.matmul(np.transpose(V),np.matmul(A,V))
    return d, V, D

