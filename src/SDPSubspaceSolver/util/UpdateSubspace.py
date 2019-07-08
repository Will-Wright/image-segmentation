import numpy as np
from scipy.linalg import orth

def UpdateSubspace(V_sub, V_new = []):
    if len(V_new) == 0:
        V_out = orth(V_sub)
    else:
        V_out = orth(np.concatenate((V_new, V_sub), axis=1))
    return V_out