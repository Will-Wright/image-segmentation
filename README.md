Description
-----------
This software package separates 

Note: this package is still in development.  

 
* solve the resulting semidefinite program with efficient software.

Contents
--------

Main function:

* `./ImSeg.py`: takes an image and an integer `k` as inputs. Returns image segmented into `k` disjoint subimages 
   - Status: not done

Source files:

* `./src/GetAdjMat.py`: computes the pixel adjacency matrix in `O(n)` flops, where `n` is the number of pixels.  The current method in `scikit-image` requires `O(n^2)`.
   - Status: not done.  Need to finish coding binary ops on adjacency relations.

* `./src/SolveSDPwithCvxopt.py`: solves semidefinite program with `cvxopt` optimization package.
   - Status: done

* `./src/SDPSubspaceSolver.py`: solves normalized cuts eigenvalue problem as a subspace SDP.  Requires 50% to 80% fewer flops than first-order methods.  
   - Status: not done

* `./test/`: contains test images and prototyping files for completing this package.


Python Dependencies
-------------------
* `Tkinter`
* `PIL`
* `cvxopt`
* `numpy`, `scipy`

