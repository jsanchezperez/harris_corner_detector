#ifndef HARRIS
#define HARRIS

#include<vector>


void harris(
  float *I,
  std::vector<int> &x,
  std::vector<int> &y,
  float alpha,
  float sigma,
  int   radius,
  float percentage,
  int   nobel_measure,
  int   nx,
  int   ny,
  int   verbose
);


#endif