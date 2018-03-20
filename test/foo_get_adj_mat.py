import Tkinter
from PIL import Image, ImageTk
import time
from sys import argv
import numpy as np
import scipy as sp
from scipy import sparse
from math import sqrt

def main(*argv):
    window = Tkinter.Tk(className="Original image")

    im = Image.open(argv[0] if len(argv) >=1 else 'mellon_collie_and_the_infinite_sadness_test.jpg').convert('LA')
    im_width = im.size[0]
    im_height = im.size[1]
    num_pixels = im_height*im_width
    max_dist_in_aff = 5

    canvas = Tkinter.Canvas(window, width=im_width, height=im_height)
    canvas.pack()

    im_arr = np.array(im).astype(np.float16)[0:im_height, 0:im_width, 0]

    im_tk = ImageTk.PhotoImage(im)
    canvas.create_image(im_width//2, im_height//2, image=im_tk)

    def callback(event):
	print "clicked at: ", event.x, event.y

    canvas.bind("<Button-1>", callback)
    Tkinter.mainloop()

    aff_arr = sp.sparse.eye(num_pixels)
    for row_shift in range(0, max_dist_in_aff-1):
        for col_shift in range(0, max_dist_in_aff-1):
            pixel_dist = sqrt(row_shift**2+col_shift**2)
            if (row_shift+col_shift) > 0 and (pixel_dist <= max_dist_in_aff):
                aff_arr = aff_arr + get_
# TODO: Resume coding affinity matrix here.  Gotta shift bands correctly.


    im_diff = get_pixel_diffs(im_arr, 1, 0)

    return (im_arr, im_diff)

def get_similarity_array(im, row_shift, col_shift):


def get_pixel_diffs(im, row_shift, col_shift):
    [r, c] = im.shape
    im_diff = np.zeros((r, c))
    im_diff[0:r-row_shift, 0:c-col_shift] = im[0:r-row_shift, 0:c-col_shift] - im[row_shift:r, col_shift:c]
    return(im_diff)

