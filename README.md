Description
-----------
This software package separates an image into `k` disjoint subimages with similar properties (color, texture, etc.).

Note: this package is still in development. When completed, it will perform the following functions:

1. Unsupervised learning of image segments


![baseball_original](baseball_original.png)

![baseball_segmented](baseball_segmented.png)


2. Semi-supervised learning of image segments

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

