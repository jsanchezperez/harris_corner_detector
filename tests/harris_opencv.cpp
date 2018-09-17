// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2018, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// All rights reserved.


#include "../src/harris.h"
#include "../src/interpolation.h"

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h> 


extern "C"
{
#include "../src/iio.h"
}


using namespace cv;
using namespace std;


/**
  *
  * Main function for computing Harris corners with OpenCV
  *
**/
  
void opencv_harris(
  float *I,        // input image
  vector<harris_corner> &corners, // output selected corners
  float k,         // Harris constant for the ....function
  float sigma_d,   // standard deviation for smoothing (image denoising)    
  float sigma_i,   // standard deviation for smoothing (pixel neighbourhood)
  float Th,        // threshold for eliminating low values
  int   strategy,  // strategy for the output corners
  int   cells,     // number of regions in the image for distributed output
  int   Nselect,   // number of output corners
  int   precision, // type of subpixel precision approximation
  int   nx,        // number of columns of the image
  int   ny,        // number of rows of the image
  int   verbose    // activate verbose mode
)
{
  int size=nx*ny;
  float *Mc = new float[size]();
  
  Mat src_gray(ny, nx, CV_32FC1, I);
  Mat dst(ny, nx, CV_32FC1 , Mc);
  
  /// Detector parameters
  int blockSize = 2;
  int apertureSize = 3;

  struct timeval start, end;
  
  if (verbose) 
  {  
     gettimeofday(&end, NULL);
     printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);

    printf(" 1-4.Computing Harris function: \t");
    gettimeofday(&start, NULL);     
  }    
  
  /// Detecting corners
  cornerHarris( src_gray, dst, blockSize, apertureSize, k, BORDER_DEFAULT );
  
  //THIS PART BELOW IS COMMOM TO OUR IMPLEMENTATION
  
  //threshold the discriminant function
  if (verbose) 
  {
     gettimeofday(&end, NULL);
     printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
     
     printf(" 5.Non-maximum suppression:  \t\t");
     gettimeofday(&start, NULL);     
  }

  //apply non-maximum suppression to select salient corners
  int radius=2*sigma_i+0.5;
  non_maximum_suppression(Mc, corners, Th, radius, nx, ny);

  if (verbose) 
  {
     gettimeofday(&end, NULL);
     printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
  }

  if (verbose) 
  {
     printf(" 6.Selecting output corners:  \t\t");
     gettimeofday(&start, NULL);     
  }

  //select output corners depending on the strategy
  //all corners; all corners sorted; N corners; distributed corners
  select_output_corners(corners, strategy, cells, Nselect, nx, ny);

  if (verbose) 
  {
     gettimeofday(&end, NULL);
     printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
  }

  if(precision==QUADRATIC_APPROXIMATION || precision==QUARTIC_INTERPOLATION)
  {
    if (verbose)
    {
       printf(" 7.Computing subpixel precision: \t");
       gettimeofday(&start, NULL);
    }

    //calculate subpixel precision
    compute_subpixel_precision(Mc, corners, nx, precision);

    if (verbose)
    {
       gettimeofday(&end, NULL);
       printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
    }
  }

  if(verbose)
    printf(" * Number of corners detected: %ld\n", corners.size());
  
  delete []Mc;
}




