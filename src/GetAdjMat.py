import Tkinter
from PIL import Image, ImageTk
import time
from sys import argv
import numpy as np
import scipy as sp
from scipy import sparse
from math import sqrt

def main(*argv):

    # Sets default values
    feat_type = 'grayscale_intensity'
    sig_feat = 1
    sig_dist = 1
    max_dist_in_aff = 1

    im = Image.open(argv[0] if len(argv) >=1 else '../test/mellon_collie_and_the_infinite_sadness_test.jpg').convert('LA')
    im_width = im.size[0]
    im_height = im.size[1]
    num_pixels = im_height*im_width

    im_arr = np.array(im).astype(np.float16)[0:im_height, 0:im_width, 0]

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
                aff_arr = aff_arr + get_similarity_bands(im_arr, row_shift, col_shift, feat_type, sig_feat, sig_dist)

    return (im_arr, aff_arr)


def get_similarity_bands(im, row_shift, col_shift, feat_type, sig_feat, sig_dist):
    rs = row_shift
    cs = col_shift
    [r, c] = im.shape
    n = r*c

    pixel_dist = sqrt(row_shift**2+col_shift**2)
    dist_exponent = -(pixel_dist/sig_dist)**2

    if feat_type == 'grayscale_intensity':
        W_feat = np.zeros((r, c))
        if (rs == 0) and (cs == 0):
            W_feat = np.ones((r, c))
        elif row_shift >= 0:
            W_feat[0:r-rs, 0:c-cs] = np.exp(-(im[0:r-rs, 0:c-cs] - im[rs:r, cs:c])**2 / sig_feat**2)
        elif row_shift < 0:
            W_feat[-rs:r, 0:c-cs] = np.exp(-(im[-rs:r, 0:c-cs] - im[0:r+rs, cs:c])**2 / sig_feat**2)

    W = np.exp(dist_exponent)*W_feat

    # Shift value will always be non-negative based on 'for' loop range
    # in main function
    shift_val = rs + r*cs
    band_arr = np.reshape(W, n, order='F')

    W_banded = sp.sparse.spdiags(band_arr, -shift_val, n, n)
    if (rs != 0) or (cs != 0):
        W_banded = W_banded + W_banded.T

    return(W_banded)




