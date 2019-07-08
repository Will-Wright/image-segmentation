# This function generates the adjacency matrix for a given image
#
# This function vectorizes the inputted image
# by mapping columns a_1, ..., a_colmax
# to v = [ [a_1], ..., [a_colmax] ] and then
# generating the adjacency matrix
# 
#                                          { exp( - ||Xi-Xj||^2 / sig_dist^2)
# W_ij = exp( -||Fi-Fj||^2 / sig_feat^2) * {         if ||Xi-Xj|| <= max_dist
#                                          { 0       else
# where Fi = feature msmt for pixel i (i.e., color scale value)
# and Xi = pixel location
# 
# The returned matrix is sparse
#


#import Tkinter
#from PIL import Image, ImageTk
#import time
from sys import argv
import numpy as np
import scipy as sp
from scipy import sparse
from math import sqrt
import warnings

def GetAdjMat(im_arr, sig_feat, sig_dist, max_dist_in_aff, is_lifted):

    # Sets default values
    feat_type = 'grayscale_intensity'
#    sig_feat = 1
#    sig_dist = 1
#    max_dist_in_aff = 1
    im_height = im_arr.shape[0]
    im_width = im_arr.shape[1]
    num_pixels = im_height*im_width
    if (im_arr.ndim == 3) and (im_arr.shape[2] == 3):
        is_RGB_image = True
    else:
        is_RGB_image = False
        
#    print("is rgb:", is_RGB_image)
#    print("im_arr.shape", im_arr.shape)
    
    if is_lifted:
        aff_arr = sp.sparse.csr_matrix((num_pixels+1, num_pixels+1))
    else:
        aff_arr = sp.sparse.csr_matrix((num_pixels, num_pixels))

    row_col_max = max_dist_in_aff+1
    for col_shift in range(0, row_col_max):
        # Prevents loop from summing duplicate similarity bands
        if col_shift == 0:
            row_min = 0
        elif col_shift > 0:
            row_min = -max_dist_in_aff
        for row_shift in range(row_min, row_col_max):
            pixel_dist = sqrt(row_shift**2+col_shift**2)
            if pixel_dist <= max_dist_in_aff:
                aff_arr = aff_arr + get_similarity_bands(im_arr, row_shift, \
                                                         col_shift, feat_type, \
                                                         sig_feat, sig_dist, \
                                                         is_lifted, is_RGB_image)

    return aff_arr


def get_similarity_bands(im_arr, row_shift, col_shift, \
                         feat_type, sig_feat, sig_dist, \
                         is_lifted, is_RGB_image):
    rs = row_shift
    cs = col_shift
    r = im_arr.shape[0]
    c = im_arr.shape[1]
    if is_RGB_image:
        dim = 3
    else:
        dim = 1
    num_pixels = r*c

    # Suppresses potential warning for `true_divide`
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Computes feature value for each pair of pixels
        if feat_type == 'grayscale_intensity':
            W_feat = np.zeros((r, c))
            if (rs == 0) and (cs == 0):
                W_feat = np.ones((r, c))
            elif row_shift >= 0:
                if is_RGB_image:
                    
                    W_feat_temp = np.zeros((r, c, dim))
                    
                    W_feat_temp[0:r-rs, 0:c-cs, :] \
                        = np.exp(  -( im_arr[0:r-rs, 0:c-cs, :] - im_arr[rs:r, cs:c, :] )**2 \
                                 / sig_feat**2)
                    
                    W_feat[0:r-rs, 0:c-cs] = np.mean(W_feat_temp[0:r-rs, 0:c-cs, :], axis=2)
                    
                    
                else:
                    W_feat[0:r-rs, 0:c-cs] \
                        = np.exp(-(im_arr[0:r-rs, 0:c-cs] - im_arr[rs:r, cs:c])**2 / sig_feat**2)
            elif row_shift < 0:
                if is_RGB_image:
                    
                    W_feat_temp = np.zeros((r, c, dim))
                    
                    W_feat_temp[-rs:r, 0:c-cs, :] \
                        = np.exp(  -( im_arr[-rs:r, 0:c-cs, :] - im_arr[0:r+rs, cs:c, :] )**2 \
                                 / sig_feat**2)
                    W_feat[-rs:r, 0:c-cs] = np.mean(W_feat_temp[-rs:r, 0:c-cs], axis=2)
                    
                    
                else:
                    W_feat[-rs:r, 0:c-cs] \
                        = np.exp(-(im_arr[-rs:r, 0:c-cs] - im_arr[0:r+rs, cs:c])**2 / sig_feat**2)

        # Scales matrix band by distance value (fixed value for each band)
        pixel_dist = sqrt(row_shift**2+col_shift**2)
        dist_exponent = -(pixel_dist/sig_dist)**2
        W = np.exp(dist_exponent)*W_feat

    # Shift value will always be non-negative based on 'for' loop range
    # in main function
    shift_val = rs + r*cs
    band_arr = np.reshape(W, num_pixels, order='F')

    # Returns symmetric matrix for all offdiagonal bands
    if is_lifted:
        W_banded = sp.sparse.spdiags(np.concatenate((band_arr, [0])), \
                                     -shift_val, num_pixels+1, num_pixels+1)
    else:
        W_banded = sp.sparse.spdiags(band_arr, -shift_val, \
                                     num_pixels, num_pixels)
    if (rs != 0) or (cs != 0):
        W_banded = W_banded + W_banded.T

    return(W_banded)




