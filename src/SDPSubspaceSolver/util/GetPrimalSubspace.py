from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import lobpcg

def GetPrimalSubspace(A, B=None, num_eigs=1, sigma = None, eig_solver = "lobpcg", \
                      eigsh_tol=1e-6, eigsh_maxiters=2000, \
                      lobpcg_tol=1e-6, lobpcg_maxiters=10000, \
                      v_prev = [], which = 'SA'):
    if eig_solver == "eigsh":
        [d, V] = eigsh(A, M=B, k = num_eigs, sigma = sigma, v0 = v_prev, \
                       which = which, tol = eigsh_tol, maxiter=eigsh_maxiters) 
    else:
        if which == 'SA':
            [d, V] = lobpcg(A, v_prev, B=B, M=None, Y=None, \
                            tol=lobpcg_tol, maxiter=lobpcg_maxiters, \
                            largest = False, verbosityLevel=0, \
                            retLambdaHistory=False, retResidualNormsHistory=False)
        else:
            [d, V] = lobpcg(-A, v_prev, B=B, M=None, Y=None, \
                            tol=lobpcg_tol, maxiter=lobpcg_maxiters, \
                            largest = False, verbosityLevel=0, \
                            retLambdaHistory=False, retResidualNormsHistory=False)

            d = -d
    return d, V
