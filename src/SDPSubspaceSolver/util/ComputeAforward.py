import numpy as np

def ComputeAforward(A,X):
    if isinstance(A,tuple):
        m = len(A)
    else:
        m = A.shape[2]

    AX = np.zeros((m,1))

    for i in range(0,m):
        if isinstance(A,tuple):
            if type(A[i]) is np.ndarray:
                AX[i,0] = np.tensordot(A[i],X,2)
            else:
                AX[i,0] = np.tensordot(A[i],X,2)
        else:
            if m == 1:
                AX[i,0] = np.tensordot(A[:,:,i],X,2)
            else:
                AX[i,0] = np.tensordot(np.squeeze(A[:,:,i]),X,2)

    return AX
