# This function solves semidefinite programs (SDPs) in the form
# with primal
#	    min   tr(C'*X)
# (SDP-P)    X
#            st   A(X) = b, X >= 0
#
# and dual
#           max   b'*y
# (SDP-D)    y
#            st   Z = C - A'(y) >= 0
#
# Note C must be a symmetric or Hermitian matrix,
# and A is a linear operator comprised of symmetric
# or Hermitian matrices Ai.
#
# This function is essentially a wrapper for cvxopt
# http://cvxopt.org/userguide/coneprog.html#semidefinite-programming
# Cvxopt has a fairly different convention for framing SDPs.
# Below is a conversion table for the SDP coefficients and variables
# for our SDP primal and dual vs the cvxopt SDP primal and dual.
#
#     our format        cvxopt format
#     ----------------------------------------------------
#           A           Gs1 = [vec(A1), ..., vec(Am)]
#           b           -c
#           C           hs1
#           X           zs1
#           y           x
#           Z           s1
#           primal      -dual
#           dual        -primal
#

from cvxopt import matrix, solvers
from numpy import array

import numpy as np
import scipy as sp

def SolveSDPwithCvxopt(A, b, C, verbose=False, gap_tol=1e-7, \
                       rel_tol=1e-6, feas_tol=1e-7):
    c = matrix(-b)
    # Note: A must be a numpy array of size (r, c, t) = (n, n, m)
    # where r = rows, c = cols, and t = no. matrices
    m = A.shape[2]
    n = A.shape[0]
    #G = [ matrix( np.transpose(np.reshape(A, (m, n*n))) ) ]
    G = [ matrix( np.reshape(A, (n*n, m)) ) ]
    h = [ matrix(C) ]
    solvers.options['show_progress'] = verbose
    solvers.options['abstol'] = gap_tol
    solvers.options['reltol'] = rel_tol
    solvers.options['feastol'] = feas_tol
    sol = solvers.sdp(c, Gs=G, hs=h)
#    print("\ny = \n")
#    print(sol['x'])
#    print("X = \n")
#    print(sol['zs'][0])
#    print("Z =\n")
#    print(sol['s1'][0])
    X = array(sol['zs'][0])
    y = array(sol['x'])
    return X, y



