#include "gaussian.h"


#include<math.h>

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
  int precision //defines the size of the window
)
{
  int i, j, k;
  float den=2*sigma*sigma;
  int size=(int) (precision*sigma)+1;
  int bdx=xdim+size;
  int bdy=ydim+size;

  //compute the coefficients of the 1D convolution kernel
  float *B=new float[size];
  for (int i=0; i<size; i++)
    B[i]=1/(sigma*sqrt(2.0*3.1415926))*exp(-i*i/den);

  float norm=0;

  //normalize the 1D convolution kernel
  for (int i=0; i<size; i++)
    norm+=B[i];

  norm*=2;
  norm-=B[0];

  for (int i=0; i<size; i++)
    B[i]/=norm;
    
  //convolution of each line of the input image
  #pragma omp parallel for
  for (k=0; k<ydim; k++)
    {
      float *R = new float[size + xdim + size]; 
      
      for (i = size; i < bdx; i++) 
        R[i] = I[k * xdim + i - size];
      
      //Dirichlet boundary conditions
      for (i = 0, j = bdx; i < size; i++, j++)
	  R[i] = R[j] = 0;
          
      for (i = size; i < bdx; i++)
        {
          float sum = B[0] * R[i];

          for (int j = 1; j < size; j++)
            sum += B[j] * (R[i - j] + R[i + j]);

          I[k * xdim + i - size] = sum;        
        }
      delete[]R;
    }

  //convolution of each column of the input image
  #pragma omp parallel for
  for (k = 0; k < xdim; k++)
    {
      float *T = new float[size + ydim + size];
      
      for (i = size; i < bdy; i++)
        T[i] = I[(i - size) * xdim + k];

      // Dirichlet boundary conditions
      for (i = 0, j = bdy; i < size; i++, j++)
	T[i] = T[j] = 0;
      
      for (i = size; i < bdy; i++)
        {
          float sum = B[0] * T[i];

          for (j = 1; j < size; j++)
            sum += B[j] * (T[i - j] + T[i + j]);

          I[(i - size) * xdim + k] = sum;
        }
        
      delete[]T;
    }
  
  delete[]B;
}
