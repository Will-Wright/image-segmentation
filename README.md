Description
-----------
This software package separates an image into `k` disjoint subimages with similar properties (color, texture, etc.).

Note: this package is still in development. When completed, it will perform the following functions:

1. Unsupervised learning of image segments: given the input image

<p align="center"> 
<img src="baseball_original.png">
<span>Here's the overlay text</span>
</p>

<p align="center">
<img src="baseball_segmented.png">
</p>

[Jianbo Shi, Jitendra Malik:
Normalized Cuts and Image Segmentation. IEEE Trans. Pattern Anal. Mach. Intell. 22(8): 888-905 (2000)](https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/papers/sm_pami00.pdf)


2. Semi-supervised learning of image segments

<p align="center">
<img src="flower_segmentation_with_contraints.png">
</p>


Contents
--------

Main function:

* `./ImSeg.py`: takes an image and an integer `k` as inputs. Segments the image into `k` disjoint subimages. Returns an array of the `k` disjoint subimages. 
   - Status: not done

Source files:

* `./src/GetAdjMat.py`: computes the pixel adjacency matrix in `O(n)` flops, where `n` is the number of pixels.  The current method in `scikit-image` requires `O(n^2)`.
   - Status: not done.  Need to finish coding binary ops on adjacency relations.

* `./src/SolveSDPwithCvxopt.py`: transforms SDP (semidefinite program) into  appropriate model format and solves it with `cvxopt` optimization package.
   - Status: done

* `./src/SDPSubspaceSolver.py`: solves normalized cuts eigenvalue problem as a subspace SDP.  Requires 50% to 80% fewer flops than first-order methods.  
   - Status: not done

* `./test/`: contains test images and prototyping files for use in finishing package.


Python Dependencies
-------------------
* `Tkinter`
* `PIL`
* `cvxopt`
* `numpy`, `scipy`

