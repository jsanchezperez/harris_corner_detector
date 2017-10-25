#ifndef INTERPOLATION
#define INTERPOLATION


/**
  *
  * Apply Newton method to find maximum of the interpolation function
  *
**/
bool maximum_interpolation(
  float *M, //values of the surfare (9 values)
  float &x, //x solution
  float &y, //y solution
  float TOL=10E-10  //stopping criterion threshold
);


#endif