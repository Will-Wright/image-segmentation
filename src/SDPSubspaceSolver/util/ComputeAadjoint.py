import numpy as np
import scipy as sp
from scipy import sparse

def ComputeAadjoint(A,y):
    m = len(A)
    n = A[0].shape[0]
    Aty = sp.sparse.csr_matrix((n,n))
    for i in range(0, m):
        if type(A[i]) is np.ndarray:
            Aty = Aty + y[i]*A[i]            
        else:                
            Aty = Aty + A[i].multiply(y[i])
    
    return Aty
#    return np.squeeze(np.dot(A,y))