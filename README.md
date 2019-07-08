Description
-----------
This software package separates an image into disjoint subimages with similar properties (color, texture, etc.). The user may also include must-link constraints.

#### New Contributions:

* more efficient method for image segmentation eigenvalue problem

* more efficient method for computing adjacency matrix


Demo Tutorial
-------------

The following results demonstrate the efficiency of our `SDPSubspaceSolver` over the previously developed `NewtonEigSolver` for image segmentation as proposed by [Anders P. Eriksson, Carl Olsson, Fredrik Kahl:
Normalized Cuts Revisited: A Reformulation for Segmentation with Linear Grouping Constraints. Journal of Mathematical Imaging and Vision 39(1): 45-61 (2011)](http://www2.maths.lth.se/vision/publdb/reports/pdf/eriksson-olsson-etal-jmiv-10.pdf)

To recreate this demo, call the function `RunDemo.main()` (dependencies listed below). 

<p align="center"> 
<img src="person_walking_small.jpg">
</p>
<p align="center">
Original image
</p>

<p align="center">
<img src="person_walking_seg.png">
</p>
<p align="center">
Segmented regions
</p>


<p align="center">
<img src="demo_fig1.png">
</p>
<p align="center">
</p>





Contents
--------

Main function:

* `./ImSeg.py`: takes an image and possible link constriants as inputs. Segments the image into disjoint subimages. Returns an array of disjoint subimages.

Key files and folders:

* `./src/SDPSubspaceSolver/*`: our new method which solves the image segmentation eigenvalue problem as a subspace SDP (semidefinite program).  Requires 50% to 80% fewer flops than the previous method by Eriksson, et al.

* `./src/NewtonEigSolver/*`: an implementation of the original Newton method for the image segmentation problem by Eriksson, et al.

* `./src/GetAdjMat.py`: computes the pixel adjacency matrix in `O(n)` flops, where `n` is the number of pixels.  The current method in `scikit-image` requires `O(n^2)`.

* `./test/`: contains test images and prototyping files.


Python Dependencies
-------------------
* `Tkinter`
* `PIL`
* `cvxopt`
* `matplotlib`, `numpy`, `scipy`


Additional Comments
-------------------

* The original non-constrained image segmentation method is [Jianbo Shi, Jitendra Malik:
Normalized Cuts and Image Segmentation. IEEE Trans. Pattern Anal. Mach. Intell. 22(8): 888-905 (2000)](https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/papers/sm_pami00.pdf)
