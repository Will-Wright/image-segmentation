# This function separates an image into meaningful, disjoint
# subimages.  The function takes as inputs a user-provided image
# and (optionally) user-selected pixels and returns
#
#
# This solver solves the normalized cuts (NCuts) model in its homogeneous, nullspace form
#
#                     x'(C + t*A2)x
# (NCuts-HN)  max min --------------
#              t   x      x'A1x
#
# The NCuts model is then solved in the matrix-lifted SDP form
#
#          min <C, X>
#           X
# (SDP-P)   st <A1, X> = 1
#              <A2, X> = 0
#
#          max  s
# (SDP-D)  s,t
#          st  C - s*A1 + t*A2 >= 0
#
# which can be solved efficiently using an iterative subspace
# method which achieves local super-quadratic convergence
#
# References:
# - for normalized cuts with grouping constraints
# --- http://www2.maths.lth.se/vision/publdb/reports/pdf/eriksson-olsson-etal-jmiv-10.pdf
# - for the subspace method used to solve (NCuts-SDP)
# --- http://homepage.divms.uiowa.edu/~oliveira/PAPER2/oliveira-stew-soma-LNCS2329-02.pdf
# - for local super-quadratic convergence rate of the subspace method
# --- https://epubs.siam.org/doi/10.1137/17M1127545
#

import tkinter
import copy
from PIL import Image, ImageTk
#import matplotlib.pyplot as plt
import scipy as sp
from sys import argv
import numpy as np
import os, sys, inspect
import time

sys.path.append(os.getcwd() + "/src")
sys.path.append(os.getcwd() + "/src/SDPSubspaceSolver")
sys.path.append(os.getcwd() + "/src/SDPSubspaceSolver/util")
sys.path.append(os.getcwd() + "/src/NewtonEigSolver")

from GetAdjMat import GetAdjMat
from SolveNCuts import SolveNCuts
from GetConstraints import GetConstraints
from GetConstraintNullspace import GetConstraintNullspace
from GetSDPModel import GetSDPModel
from SDPSubspaceSolver import SDPSubspaceSolver
from NewtonEigSolver import NewtonEigSolver
from GetPrimalSubspace import GetPrimalSubspace


def ImSeg(im_path = 'test/amoeba.jpg', Tk_im_file = [], \
          attempt_resize=True, resize_ratio=0.5, \
          resize_min_num_pix=200, \
          sig_dist=5, sig_feat=2, max_dist=5, \
          pr_tol=1e-8, du_tol=1e-8, gap_tol=1e-6, \
          SDPSS_max_iters=10, SDPSS_pr_update=1, eig_solver = "lobpcg", \
          eigsh_tol_init=1e-3, eigsh_tol_update=1e-1, eigsh_tol_min=1e-6, \
          lobpcg_tol_init=1e-4, lobpcg_tol_update=1e-1, lobpcg_tol_min=1e-9, \
          cvx_verbose=False, must_link_list1=[], must_link_list2=[], \
          use_RGB = True, \
          use_NewtonEigSolver=False, Newton_max_iters=10, \
          Newton_eig_solver="eigsh", use_Newton_step=True, \
          is_demo_mode=False, return_data=False, \
          verbosity=1):

    # Sets parameters for resized problems
    resize_filter = Image.ANTIALIAS
#    resize_filter = Image.NEAREST
#    resize_filter = Image.BILINEAR
#    resize_filter = Image.BICUBIC

    # Opens image
    if not Tk_im_file:
        im_RGB = Image.open(im_path)
    else:
        im_RGB = Tk_im_file

    im_width = im_RGB.size[0]
    im_height = im_RGB.size[1]
    im_num_pix = im_width*im_height

    # Gets linking constraints for image
    (must_link_list1, must_link_list2) \
        = GetConstraints(im_width, im_height, im_RGB, \
                         must_link_list1, must_link_list2, \
                         is_demo_mode)

    # Solves im seg subproblem(s) to initialize full problem
    if attempt_resize:

        resize_dims_list = list()
        resize_ratio_list = list()
        resize_dims_list.append((im_width, im_height))
        resize_ratio_list.append(1)
        resize_ratio_list.append(resize_ratio)
        im_num_pix_temp = resize_dims_list[-1][0]*resize_dims_list[-1][1]

        while (im_num_pix_temp >= resize_min_num_pix):
            im_size_resize = (int(round(resize_ratio*resize_dims_list[-1][0])), \
                              int(round(resize_ratio*resize_dims_list[-1][1])))
            resize_dims_list.append(im_size_resize)
            resize_ratio_list.append(resize_ratio**(len(resize_dims_list)))
            im_num_pix_temp = im_num_pix * (resize_ratio_list[-1]**2)
        resize_ratio_list.pop()

        while len(resize_ratio_list) > 1:
            im_size_resize = resize_dims_list.pop()
            resize_ratio = resize_ratio_list.pop()
            im_RGB_resize = im_RGB.resize(im_size_resize, resize_filter)
            if verbosity>=1:
                print("Image resized to", im_size_resize)

            if (len(must_link_list1) + len(must_link_list2) >= 1):
                mll_resize1 = [[int(resize_ratio*i) for i in tups] \
                             for tups in must_link_list1]
                mll_resize2 = [[int(resize_ratio*i) for i in tups] \
                             for tups in must_link_list2]
            else:
                mll_resize1 = []
                mll_resize2 = []

            [im_RGB_arr, adj_mat, must_link_list1, must_link_list2, \
             X_sub, V_sub, im, im_RGB, N, data] = ImSeg(\
                        im_path = [], Tk_im_file = im_RGB_resize, \
                        attempt_resize=False, \
                        resize_min_num_pix=resize_min_num_pix, \
                        cvx_verbose=cvx_verbose, eig_solver = eig_solver, \
                        pr_tol=pr_tol, du_tol=du_tol, gap_tol=gap_tol, \
                        sig_dist=sig_dist, sig_feat=sig_feat, max_dist=max_dist, \
                        eigsh_tol_init=eigsh_tol_init, \
                        eigsh_tol_update=eigsh_tol_update, \
                        eigsh_tol_min=eigsh_tol_min, \
                        lobpcg_tol_init=lobpcg_tol_init, \
                        lobpcg_tol_update=lobpcg_tol_update, \
                        lobpcg_tol_min=lobpcg_tol_min, \
                        SDPSS_max_iters=SDPSS_max_iters, \
                        SDPSS_pr_update=SDPSS_pr_update, \
                        must_link_list1 = mll_resize1,  \
                        must_link_list2 = mll_resize2, \
                        use_RGB=use_RGB, \
                        use_NewtonEigSolver=use_NewtonEigSolver, \
                        Newton_eig_solver=Newton_eig_solver, \
                        use_Newton_step=use_Newton_step, \
                        Newton_max_iters=Newton_max_iters, \
                        is_demo_mode=is_demo_mode, return_data=True, \
                        verbosity=verbosity)[6];

    if (len(must_link_list1) + len(must_link_list2) >= 1):
        is_lifted = True
        if verbosity>=1:
            print("Image segmentation (NCuts) called with link constraints")
            print("   for image of size", im_RGB.size)
        start = time.process_time()
        N = GetConstraintNullspace(im_height, im_width, \
                                   must_link_list1, must_link_list2)
        end = time.process_time()
        GCN_time = float(end) - float(start)
        if verbosity>=1:
            print("GetConstraintNullspace time  : ", "%4.2f" % GCN_time)

    else:
        is_lifted = False
        if verbosity>=1:
            print("Image segmentation (NCuts) called with no link constraints")
        N = []

    start = time.process_time()
    im_RGB_arr = np.array(im_RGB).astype(np.float32)
    adj_mat = GetAdjMat(im_RGB_arr, sig_feat, sig_dist, max_dist, is_lifted)
    end = time.process_time()
    GAM_time = float(end) - float(start)
    if verbosity>=1:
        print("GetAdjMat time               : ", "%4.2f" % GAM_time)

    start = time.process_time()
    if is_lifted:
        (A, b, C, V_init) = GetSDPModel(adj_mat, N)
        SDPSS_out = SDPSubspaceSolver(A, b, C, V_init=V_init, \
                 pr_tol=pr_tol, du_tol=du_tol, gap_tol=gap_tol, \
                 eig_solver = eig_solver, \
                 eigsh_tol_init=eigsh_tol_init, \
                 eigsh_tol_update=eigsh_tol_update, \
                 eigsh_tol_min=eigsh_tol_min, \
                 eigsh_maxiters = 100000, \
                 lobpcg_tol_init=lobpcg_tol_init, \
                 lobpcg_tol_update=lobpcg_tol_update, \
                 lobpcg_tol_min=lobpcg_tol_min, \
                 lobpcg_maxiters = 100000, \
                 max_iters=SDPSS_max_iters, \
                 cvx_verbose=cvx_verbose, \
                 pr_dim_update=SDPSS_pr_update, \
                 return_data=return_data, \
                 verbosity=verbosity)

        if return_data:
            (X_sub, V_sub, y, SDPSS_data) = SDPSS_out
        else:
            (X_sub, V_sub, y) = SDPSS_out
        [d, V] = np.linalg.eig(X_sub)
        v_SDPSS = N.dot(V_sub.dot(V[:,0]))[0:im_width*im_height]
        v_SDPSS = v_SDPSS/np.linalg.norm(v_SDPSS)

        if use_NewtonEigSolver:
            NES_out = NewtonEigSolver(A, b, C, V_init=V_init, \
                 pr_tol=pr_tol, du_tol=du_tol, gap_tol=gap_tol, \
                 max_iters=Newton_max_iters, \
                 eig_solver = Newton_eig_solver, \
                 use_Newton_step=use_Newton_step, \
                 eigsh_tol_init=eigsh_tol_init, \
                 eigsh_tol_update=eigsh_tol_update, \
                 eigsh_tol_min=eigsh_tol_min, \
                 eigsh_maxiters = 100000, \
                 lobpcg_tol_init=lobpcg_tol_init, \
                 lobpcg_tol_update=lobpcg_tol_update, \
                 lobpcg_tol_min=lobpcg_tol_min, \
                 lobpcg_maxiters = 100000, \
                 return_data=return_data, \
                 verbosity=verbosity)
            if return_data:
                (V1, f, NES_data) = NES_out
            else:
                (V1, f) = NES_out
            v_Newton = N.dot(V1[:,0])[0:im_width*im_height]
            v_out = v_Newton

        else:
            v_out = v_SDPSS

        end = time.process_time()

    else:
        (d, v_NCuts) = SolveNCuts(adj_mat, eigsh_tol_min)
        v_out = v_NCuts
        end = time.process_time()

    # Converts NCuts solution to image
    v_out_sign = np.sign(v_out)
    im_seg_Newton_arr = np.reshape(v_out_sign, (im_width, im_height)).T
    im_seg_RGB_arr = np.zeros((im_height, im_width, 3), 'uint8')
    for i in range(0,im_height):
        for j in range(0,im_width):
            if im_seg_Newton_arr[i,j] == 1:
                im_seg_RGB_arr[i,j, 0] = 255
                im_seg_RGB_arr[i,j, 1] = 153
                im_seg_RGB_arr[i,j, 2] = 51
            else:
                im_seg_RGB_arr[i,j, 0] = 153
                im_seg_RGB_arr[i,j, 1] = 0
                im_seg_RGB_arr[i,j, 2] = 153

    im_seg_RGB = Image.fromarray(im_seg_RGB_arr)

    if is_demo_mode:
        # Displays segmented image
        time.sleep(3)
        im_seg_RGB.show(title="Newton Image")

    if return_data:
        if is_lifted:
            if use_NewtonEigSolver:
                data = {'SDPSS_data': SDPSS_data,
                        'NES_data': NES_data}
            else:
                data = {'SDPSS_data':SDPSS_data}
        else:
            X_sub = {}
            V_sub = {}
            data = {}
        out = (im_RGB_arr, adj_mat, \
               must_link_list1, must_link_list2, \
               X_sub, V_sub, im_seg_RGB_arr, \
               im_RGB, N, data)
    else:
        out = im_seg_RGB_arr

    return out

