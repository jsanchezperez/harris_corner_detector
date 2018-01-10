#ifndef _GAUSSIAN_H
#define _GAUSSIAN_H


#define STD_GAUSSIAN 0
#define FAST_GAUSSIAN 1
#define NO_GAUSSIAN 2

/**
 *
 * Convolution with a Gaussian 
 *
 */
void gaussian(
  float *I,    //input/output image
  int   nx,    //image width
  int   ny,    //image height
  float sigma, //Gaussian sigma
  int   type=FAST_GAUSSIAN, //type of Gaussian convolution 
  int   K=3    //defines the number of iterations or window precision
);


#endif
