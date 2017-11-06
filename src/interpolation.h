#ifndef INTERPOLATION
#define INTERPOLATION


#include "harris.h"


/**
  *
  * Apply Newton method to find maximum of the interpolation function
  *
**/
bool maximum_interpolation(
  float *M, //values of the surfare (9 values)
  harris_corner &corner, //corner
  float TOL=10E-10  //stopping criterion threshold
);


#endif