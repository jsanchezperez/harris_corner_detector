#ifndef _GRADIENT_H
#define _GRADIENT_H


#define CENTRAL_DIFFERENCES 0
#define SOBEL_OPERATOR 1

/**
  *
  * Function to compute the gradient
  *
**/
void gradient(
  float *I,  //input image
  float *dx, //computed x derivative
  float *dy, //computed y derivative
  int   nx,  //image width
  int   ny,  //image height
  int   type //type of gradient
);

#endif
