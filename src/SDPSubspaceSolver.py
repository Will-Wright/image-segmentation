# This function solves semidefinite programs (SDPs)
# using a subspace method on the primal variable X.
#
# This function uses the following SDP form with primal
#           min   tr(C'*X)
# (SDP-P)    X
#            st   A(X) = b, X >= 0
#
# and dual
#           max   b'*y
# (SDP-D)    y
#            st   Z = C - A'(y) >= 0
#
# A is a linear operator defined by
# 	forward: A(X) = vec( tr(A1'*x), ..., tr(Am'*X) )
# 	adjoint: A'(y) = y1*A1 + ... + ym*Am
#	Ai are all n-by-n symmetric or Hermitian
# C is a symmetric or Hermitian matrix
# b is a real m-vector
# X is the primal variable, an n-by-n symmetric or Hermitian matrix
# y is the dual variable, a real m-vector
#
# Note: this subspace method is particularly efficient for
# SDPs with a very large objective matrix (e.g., n > 1000)
# and very few primal constraints (e.g., m < 50).
# Well-suited problems include trace-ratio problems, image segmentation
# with normalized cuts, and trust-region subproblems.
#
#
# The following papers are notable in developing this subspace
# method and its theory.
#
# 2002 Oliveira Stewart Soma, A Subspace Semidefinite Programming
# for Spectral Graph Partitioning
# 	- discusses this subspace method for first time in literature
#	- presents a few small applications, no theory
#
# 2015 Kangal Meerbergen Mengi Michiels, A Subspace Method
# for Large Scale Eigenvalue Optimization
#	- reintroduces subspace SDP method as eigenvalue method
#	- establishes basic convergence theory for general case
#
# 2017 Kressner Lu Vandereycken, Subspace acceleration
# for the Crawford number and related eigenvalue optimization problems
#	- establishes super-quadratic (1 + sqrt(2)) local convergence
#  	  for subspace method with m = 1 (only 1 primal constraint)
#

