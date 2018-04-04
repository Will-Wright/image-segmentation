import Tkinter
from PIL import Image, ImageTk
import time
from sys import argv
import numpy as np
import scipy as sp
from scipy import sparse
from math import sqrt

def main(*argv):
#    window = Tkinter.Tk(className="Original image")

    # Sets default values
    feature_type = 'grayscale_intensity'
    sigma_feature = 1
    sigma_distance = 1
    max_dist_in_aff = 3

    im = Image.open(argv[0] if len(argv) >=1 else '../test/mellon_collie_and_the_infinite_sadness_test.jpg').convert('LA')
    im_width = im.size[0]
    im_height = im.size[1]
    num_pixels = im_height*im_width

#    canvas = Tkinter.Canvas(window, width=im_width, height=im_height)
#    canvas.pack()

    im_arr = np.array(im).astype(np.float16)[0:im_height, 0:im_width, 0]

#    im_tk = ImageTk.PhotoImage(im)
#    canvas.create_image(im_width//2, im_height//2, image=im_tk)

#    def callback(event):
#	print "clicked at: ", event.x, event.y

#    canvas.bind("<Button-1>", callback)
#    Tkinter.mainloop()

    aff_arr = sp.sparse.eye(num_pixels)

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
                1+1


    return (im_arr)


def get_similarity_bands(im, row_shift, col_shift, feat_type, sig_feat, sig_dist):
    [r, c] = im.shape
    n = r*c

    pixel_dist = sqrt(row_shift**2+col_shift**2)
    dist_exponent = -(pixel_dist/sig_dist)**2

    if feat_type == 'grayscale_intensity':
    im_diff_squared = np.zeros((r, c))
        if (row_shift == 0) and (col_shift == 0):
            # Trivial case
        elif row_shift >= 0:
            im_diff_squared[0:r-row_shift, 0:c-col_shift] = (im[0:r-row_shift, 0:c-col_shift] - im[row_shift:r, col_shift:c])**2
        elif row_shift < 0:
            im_diff_squared[-row_shift:r, 0:c-col_shift] = (im[-row_shift:r, 0:c-col_shift] -         im[0:r+row_shift, col_shift:c])**2


    # TODO:
    # 1) Create matrix M = exp(-dist - affinity)
    #   Make sure appropriate entries in M remain ZEROS
    # 2) Get code to compile, run, and then write simple tests to verify
    #   Note: shouldn't need to revise (much), already checked
    #         that M RESHAPES CORRECTLY to M_banded


    # Shift value will always be non-negative based on for loop range
    # in main function
    shift_val = row_shift + r*col_shift
    band_shifts = [shift_val, -shift_val]
    band_array = np.reshape(M, n, order='F')
    bands = [band_array, band_array]

    M_banded = sp.sparse.spdiags(bands, band_shifts, n, n).toarray()

    return(M_banded)




