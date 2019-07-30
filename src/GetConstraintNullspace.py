# Computes nullspace matrix to enforce must-link, cannot-link constraints
#
# Returns matrix `N` which transforms the Laplacian matrix `M = D - W`
# and diagonal matrix `D` into smaller sparse matrices `A1 = N'*M*N`
# and `A2 = N'*D*N` which enforce link constraints
#
# The desired solution is then set to `x = N*y`

import numpy as np
import scipy as sp
from scipy import sparse
import warnings

def get_pixel_vec_idx(pix_wd_idx, pix_ht_idx, im_height):
    return pix_ht_idx + pix_wd_idx*im_height

def get_link_array(link_tuple, im_height):
    link_arr = list()
    for i in range(0, len(link_tuple)):
        (pix_wd_idx, pix_ht_idx) = link_tuple[i]
        link_arr.append(get_pixel_vec_idx(pix_wd_idx, pix_ht_idx, im_height))
    link_arr.sort()
    # Remove duplicates
    link_arr = list(dict.fromkeys(link_arr))
    return link_arr

def GetConstraintNullspace(im_height, im_width, must_link_tuple1, must_link_tuple2):

    num_pixels = im_height*im_width
    link_arr1 = get_link_array(must_link_tuple1, im_height)
    link_arr2 = get_link_array(must_link_tuple2, im_height)

    # Construct N using these arrays
    num_cons = len(link_arr1) + len(link_arr2)
    full_dim = num_pixels + 1
    N_rank = full_dim - num_cons

    # Build last column of N
    last_v_row_idx = np.concatenate((link_arr1, link_arr2, [full_dim-1]))
    last_v_col_idx = (N_rank-1)*np.ones(num_cons+1,)
    last_v = (1/np.sqrt(1+num_cons))*np.concatenate((-1*np.ones(len(link_arr1)),\
                                                     np.ones(len(link_arr2)), [1] ))

    row = np.concatenate((np.delete(np.arange(0,full_dim), \
                                    np.concatenate((link_arr1, \
                                                    link_arr2, \
                                                    [full_dim-1])) ), \
                          last_v_row_idx))
    col = np.concatenate((np.arange(0,N_rank-1), last_v_col_idx))
    data = np.concatenate((np.ones(N_rank-1,), last_v))

    N = sp.sparse.csr_matrix((data, (row, col)), shape=(full_dim, N_rank))

    return N
