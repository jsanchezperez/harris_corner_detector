#include "gradient.h"

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
)
{
    #pragma omp parallel for
    for(int i=0; i<nx; i++)
    {
      dx[i]=dy[i]=0;
      dx[nx*(ny-1)+i]=dy[nx*(ny-1)+i]=0;      
    }
    
    #pragma omp parallel for
    for(int i=1; i<ny-1; i++)
    {
        int k=i*nx;
	dx[k]=dy[k]=0;
        for(int j=1; j<nx-1; j++)
        {
	    k++;
            dx[k]=0.5*(input[k+1]-input[k-1]);
            dy[k]=0.5*(input[k+nx]-input[k-nx]);
        }
        dx[k+1]=dy[k+1]=0;
    }
}
