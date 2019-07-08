import numpy as np

def DetermineIfConverged(pr_obj, du_obj, pr_res, du_res, pr_tol, du_tol, gap_tol):
    converged = (np.absolute(pr_obj - du_obj) <= gap_tol) and (pr_res <= pr_tol) and (du_res + du_tol >= 0)
    return converged
