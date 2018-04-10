# This function separates an image into meaningful, disjoint
# subimages.  The function takes as inputs a user-provided image
# and (optionally) user-selected pixels and returns
# TODO: What to return???
#
# TODO: In README, put clear list of terminal instructions for all package installations.
#
# This solver solves the following normalized cuts (Ncuts) model
# [Eriksson, 2007]
#
# 	           x'(A1 + t*A3)x
# (Ncuts)  max min --------------
# 	    t   x      x'A2x
#
# Here, the Ncuts model is considered in an equivalent SDP form
#
# 		max  s
# (Ncuts-SDP)	s,t
#  		 st  A1 - s*A2 + t*A3 >= 0
#
# which can be solved efficiently using an iterative subspace
# method which achieves local super-quadratic convergence
# [Kressner Lu Vandereycken, 2017].
#

import Tkinter
from PIL import Image, ImageTk
from sys import argv
import numpy as np

def main(**kwargs):
    if 'image_path' in kwargs:
        im_path = kwargs['image_path']
        print 'Image path provided by user: ', im_path
    else:
        im_path = 'test/mellon_collie_and_the_infinite_sadness_test.jpg'
        print 'Default image chosen from path: ', im_path
    window = Tkinter.Tk(className='Initial image')

# TODO: Remove conversion to grayscale
    # Opens image
    im = Image.open(im_path).convert('LA') # .convert('LA') makes image grayscale
    im_width = im.size[0]
    im_height = im.size[1]
    canvas = Tkinter.Canvas(window, width=im_width, height=im_height)
    canvas.pack()
    im_tk = ImageTk.PhotoImage(im)

    # Displays image and allows user to select pixels for image subregions
    canvas.create_image(im_width//2, im_height//2, image=im_tk)
    region1_array = np.zeros((100, 2))

# TODO NEXT: Correctly store pixel coords

#    arr_idx = 0
    def callback(event):
#        global arr_idx
#        region1_array[arr_idx, :] = [event.x, event.y]
        print 'clicked at: ', event.x, event.y
#        arr_idx += 1
    canvas.bind("<Button-1>", callback)
    Tkinter.mainloop()

    # Converts image to numpy array   TODO: Remove grayscale biz
    im_array = np.array(im).astype(np.float16)[0:im_height, 0:im_width, 0]

    return region1_array



