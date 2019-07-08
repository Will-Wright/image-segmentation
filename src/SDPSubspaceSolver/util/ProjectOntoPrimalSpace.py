import numpy as np

def ProjectOntoPrimalSpace(V_sub, X_sub):
    return np.dot(np.dot(V_sub,X_sub), V_sub.transpose())