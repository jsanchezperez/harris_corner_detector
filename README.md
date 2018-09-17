Harris Corner Detector
======================

*******
SUMMARY
*******

This program implements the Harris corner detector. This feature detector 
relies on the analysis of the eigenvalues of the auto-correlation matrix.
The algorithm comprises seven steps, including several measures for the 
classification of corners, a generic non-maximum suppression method for 
selecting interest points, and the possibility to obtain corners position with 
subpixel accuracy.

This program is part of an IPOL publication:
http://www.ipol.im/

Reference articles:

[1] J. Sánchez, N. Monzón and A. Salgado, "An Analysis and Implementation of 
the Harris Corner Detector", Image Processing On line, Preprint (2018)


******
AUTHOR
******

Javier Sánchez Pérez <jsanchez@ulpgc.es> 
Centro de Tecnologías de la Imagen (CTIM) 
Universidad de Las Palmas de Gran Canaria


*******
VERSION
*******

Version 1, released on September 20, 2018


*******
LICENSE
*******

This program is free software: you can use, modify and/or redistribute it
under the terms of the simplified BSD License. You should have received a
copy of this license along this program. If not, see
<http://www.opensource.org/licenses/bsd-license.html>.

Copyright (C) 2018, Javier Sánchez Pérez <jsanchez@ulpgc.es>
All rights reserved.


***********
COMPILATION
***********

Required environment: Any unix-like system with a standard compilation
environment (make and C/C++ compilers)

Required libraries: libpng, lipjpeg, libtiff

Compilation instructions: run "make" to produce an executable


*****
USAGE
*****

The program reads an input images, take some parameters and produce a list of
interest points. The meaning of the parameters is thoroughly discussed on the 
accompanying IPOL article. Usage instructions:

  Usage: bin/harris_corner_detector image [OPTIONS] 

  Harris corner detector:
  'image' is an input image to detect features on.
  -----------------------------------------------
  OPTIONS:
  --------
   -o name  output image with detected corners 
   -f name  write points to file
   -z N     number of scales for filtering out corners
              default value 1
   -s N     choose smoothing: 
              0.precise Gaussian; 1.fast Gaussian; 2.no Gaussian
              default value 1
   -g N     choose gradient: 
              0.central differences; 1.Sobel operator
              default value 0
   -m N     choose measure: 
              0.Harris; 1.Shi-Tomasi; 2.Harmonic Mean
              default value 0
   -k N     Harris' K parameter
              default value 0.060000
   -d N     Gaussian standard deviation for derivation
              default value 1.000000
   -i N     Gaussian standard deviation for integration
              default value 2.500000
   -t N     threshold for eliminating low values
              default value 130
   -q N     strategy for selecting the output corners:
              0.all corners; 1.sort all corners;
              2.N corners; 3.distributed N corners
              default value 0
   -c N     regions for output corners (1x1, 2x2,...NxN):
              default value 3
   -n N     number of output corners
              default value 2000
   -p N     subpixel accuracy
              0.no subpixel; 1.quadratic approximation; 2.quartic interpolation
              default value 1
   -v       switch on verbose mode 



Execution examples:

  1.Default parameters and write corners to an output file:
    
   >bin/harris_corner_detector data/building.png -f corners.txt 
  
  2.Using Shi-Tomasi measure, threshold=10 and select best 1000 output corners:
    
   >bin/harris_corner_detector data/building.png -m 1 -t 10 -q 2 -n 1000 \
    -f corners.txt -o corners.png -v
   
   
If a parameter is given an invalid value it will take a default value.


*************
LIST OF FILES
*************

Directory src:
--------------
gaussian.cpp	: Different implementions of a convolution with a Gaussian function
gradient.cpp: Functions for computing the gradient of an image
harris.cpp:   This is the main program that implements the Harris method
iio.c:        Functions for reading and writing images 
interpolation.cpp: Functions for computing corners with sub-pixel accuracy
main.cpp:     Main algorithm to read parameters and write results
zoom.cpp:     Function for zooming out images by a factor of 2

Directory tests:
----------------
This directory contains the code for generating the tests included in the 
reference article.

harris_opencv.cpp: Implementation of the Harris method using the OpenCV code
bicubic_interpolation.cpp: Allows to transform an image through a parametric 
  model (scale, rotation, affinity, homography) using bicubic interpolation
test.cpp: Executes all the tests in the paper and generates data for graphics
graphics: directory for generating the graphics with gnuplot

