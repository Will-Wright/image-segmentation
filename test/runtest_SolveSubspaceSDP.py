# This function verifies the function SolveSubspaceSDP
# is correctly transforming SDPs from our form to
# the cvxopt form.
#
# TODO: We can see cvxopt/solvers.sdp correctly solves (A, b, C) form.
# Now we must code numerical verification and print "Test status: PASS"
#

import numpy as np
import sys
sys.path.insert(0, "../src")
import SolveSubspaceSDP

def main():
    A = np.zeros((3, 3, 2))
    A[:,:, 0] = np.array( [ [2., 6., 4.], [6., -2., 7.], [4., 7., 10.]  ] )
    A[:,:, 1] = np.array( [ [1., 7., 6.], [7., 14., 9.], [6., 9., 8.]  ] )
    b = np.array([[9.], [15.]])
    C = np.array( [ [5., 2., 1.], [2., 7., 3.], [1., 3., 8.] ] )
    X = np.array( [ [0.152513984269634, 0.206811764229985, 0.175877250073743],
            [0.206811764229985, 0.280440581743945, 0.238492787099002],
            [0.175877250073743, 0.238492787099002, 0.202819500417786] ] )
    y = np.array( [ [0.183267184630601], [0.353917564707084] ] )
    Z = np.array( [ [4.279548066242836, -1.577026060733195, -1.856574126764909],
            [-1.577026060733195, 2.411688463573145, -1.468128374777965],
            [-1.856574126764909, -1.468128374777965, 3.335987636248442] ] )
    SolveSubspaceSDP.main(A, b, C)


