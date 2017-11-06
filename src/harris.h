#ifndef _HARRIS_H
#define _HARRIS_H

#include<vector>

//Measures for the Harris discriminant function
#define HARRIS_MEASURE 0
#define SHI_TOMASI_MEASURE 1
#define HARMONIC_MEAN_MEASURE 2

//Strategies for selecting output corners
#define ALL_CORNERS 0
#define N_CORNERS 1
#define DISTRIBUTED_N_CORNERS 2


/**
  *
  * Structure for handling Harris' corners
  *
**/
struct harris_corner{
  float x,y; //position of the corner
  float Mc;  //Harris measure for this corner
} ;


/**
  *
  * Function for computing Harris corners
  *
**/
void harris(
  float *I,         // input image
  std::vector<harris_corner> &corners, // output selected corners
  int   measure,    // measure for the discriminant function
  float k,          // Harris constant for the ....function
  float sigma_i,    // standard deviation for smoothing (image denoising)    
  float sigma_n,    // standard deviation for smoothing (pixel neighbourhood)
  int   radius,     // radius of the autocorralation matrix
  int   select_strategy, // strategy for the output corners
  int   Nselect,    // number of output corners
  int   precision,  // enable subpixel precision
  int   nx,         // number of columns of the image
  int   ny,         // number of rows of the image
  int   verbose,    // activate verbose mode
  int   forensics   // activate forensics mode  
);

#endif
