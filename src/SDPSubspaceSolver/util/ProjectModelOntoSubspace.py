import numpy as np

def ProjectModelOntoSubspace(A, C, V_sub):
    C_sub = V_sub.transpose().dot(C.dot(V_sub))
    n_sub = V_sub.shape[1]
    m = len(A)
    A_sub = np.zeros((n_sub, n_sub, m))
    for j in range(0, m):
        A_sub[:,:,j] = V_sub.transpose().dot(A[j].dot(V_sub))
    return A_sub, C_sub