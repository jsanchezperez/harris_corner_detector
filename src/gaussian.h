#ifndef _GAUSSIAN_H
#define _GAUSSIAN_H

/**
 *
 * Convolution with a Gaussian
 *
 */
void
gaussian (
  float *I,     //input/output image
  int xdim,     //image width
  int ydim,     //image height
  float sigma,  //Gaussian sigma
  int precision=5 //defines the size of the window
);


/**
 *
 * Convolution with a Gaussian using Stacked Integral Images
 *
 */
void
gaussian_sii (
  float *I,     //input/output image
  int   nx,     //image width
  int   ny,     //image height
  float sigma,  //Gaussian sigma
  int   K=5     //defines the number of iterations
);


#endif
