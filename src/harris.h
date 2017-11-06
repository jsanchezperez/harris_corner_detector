#ifndef _HARRIS_H
#define _HARRIS_H

#include<vector>


#define HARRIS_MEASURE 0
#define SHI_TOMASI_MEASURE 1
#define HARMONIC_MEAN_MEASURE 2


void harris(
  float *I,
  std::vector<float> &x,
  std::vector<float> &y,
  float alpha,
  float sigma_i,
  float sigma_n,
  int   radius,
  float percentage,
  int   measure,
  int   precision,
  int   nx,
  int   ny,
  int   verbose,
  int   forensics
);

#endif
