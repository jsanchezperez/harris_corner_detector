#ifndef INTERPOLATION
#define INTERPOLATION

#define QUADRATIC_APPROXIMATION 0
#define QUARTIC_INTERPOLATION   1

/**
  *
  * Find maximum of a quadratic approximation function
  *
**/
bool quadratic_approximation(
  float *M,  //values of the interpolation function (9 values)
  float &x,  //corner x-position
  float &y,  //corner y-position
  float &Mo  //maximum of the interpolation function
);

/**
  *
  * Find maximum of a quartic interpolation function
  *
**/
bool quartic_interpolation(
  float *M,  //values of the interpolation function (9 values)
  float &x,  //corner x-position
  float &y,  //corner y-position
  float &Mo, //maximum of the interpolation function
  float TOL=1E-10  //stopping criterion threshold
);


#endif
