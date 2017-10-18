#ifndef _GRADIENT_H
#define _GRADIENT_H

/**
  *
  * Function to compute the gradient with centered differences
  *
**/
void gradient(
    const float *input,  //input image
    float *dx,           //computed x derivative
    float *dy,           //computed y derivative
    const int nx,        //image width
    const int ny         //image height
);

#endif
