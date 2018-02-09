# This function separates an image into meaningful, disjoint
# subimages.  The function takes as inputs a user-provided image
# and (optionally) user-selected pixels and returns
# TODO: What to return???
#
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

from cv2 import *
import numpy as np

def main(image_file):
    image = cv2.imread(image_file)
    image_gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    printBanner()





def printBanner():
    print("hi, bitch")








