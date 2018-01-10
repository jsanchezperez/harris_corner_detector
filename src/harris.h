#ifndef _HARRIS_H
#define _HARRIS_H

#include<vector>

//Measures for the Harris discriminant function
#define HARRIS_MEASURE 0
#define SHI_TOMASI_MEASURE 1
#define HARMONIC_MEAN_MEASURE 2

//Strategies for selecting output corners
#define ALL_CORNERS 0
#define ALL_CORNERS_SORTED 1
#define N_CORNERS 2
#define DISTRIBUTED_N_CORNERS 3


/**
  *
  * Structure for handling Harris' corners
  *
**/
struct harris_corner{
  float x,y; //position of the corner
  float Mc;  //Harris measure for this corner
};


/**
  *
  * Main function for computing Harris corners
  *
**/
void harris(
  float *I,        // input image
  std::vector<harris_corner> &corners, // output selected cornersç
  int   gauss,     // type of Gaussian 
  int   grad,      // type of gradient
  int   measure,   // measure for the discriminant function
  float k,         // Harris constant for the ....function
  float sigma_d,   // standard deviation for smoothing (image denoising)    
  float sigma_i,   // standard deviation for smoothing (pixel neighbourhood)
  float Th,        // threshold for eliminating low values
  int   strategy,  // strategy for the output corners
  int   cells,     // number of regions in the image for distributed output
  int   Nselect,   // number of output corners
  int   precision, // enable subpixel precision
  int   nx,        // number of columns of the image
  int   ny,        // number of rows of the image
  int   verbose    // activate verbose mode
);

#endif
