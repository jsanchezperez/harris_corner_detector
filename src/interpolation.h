#ifndef INTERPOLATION
#define INTERPOLATION

/**
  *
  * Apply Newton method to find maximum of the interpolation function
  *
**/
bool maximum_interpolation(
  float *M,  //values of the interpolation function (9 values)
  float &x,  //corner x-position
  float &y,  //corner y-position
  float &Mo, //maximum of the interpolation function
  float TOL=1E-10  //stopping criterion threshold
);

#endif