#include <stdio.h>
#include <math.h>
#include <float.h>
#include <sys/time.h> 

extern "C"
{
#include "iio.h"
}

#include "harris.h"
#include "gradient.h"
#include "gaussian.h"


// ........................

/**
  *
  * Function for computing Harris's corner measure
  *
**/
void harris_measure(
  float *A, float *B, float *C,
  float *Mc,
  int nx, int ny, 
  float k,
  float &max, float &min
)
{
  int size = nx*ny;

  for (int i=0;i<size;i++)
  {
      float detA   = A[i]*C[i] - B[i]*B[i];
      float traceA = A[i] + C[i];

      Mc[i] = detA - k*traceA*traceA;
      
      if (Mc[i] > max)  max = Mc[i];
      if (Mc[i] < min)  min = Mc[i];
  }
}



/**
  *
  * Function for computing Zseliski measure
  *
**/
void zseliski_measure(
  float *A, float *B, float *C,
  float *Mc,
  int nx, int ny,
  float &max, float &min
)
{
  int size = nx*ny;

  for (int i=0;i<size;i++)
  {
      float detA   = A[i]*C[i] - B[i]*B[i];
      float traceA = A[i] + C[i];

      Mc[i] = 2*detA / (traceA+0.0001);

      if (Mc[i]>max)  max = Mc[i];
      if (Mc[i]<min)  min = Mc[i];
  }
}



/**
  *
  * Function for computing Shid-Tomasi measure
  *
**/
void shi_tomasi_measure(
  float *A, float *B, float *C,
  float *Mc,
  int nx, int ny, 
  float th
)
{
  printf("\n\n Shidtomasi_measure: not implemented yet !!! \n\n");
}




/**
  *
  * Function for computing the spiral order indexes
  *
**/
void spiral_order(
  int *index,
  int radius,
  int nx,
  int ny
)
{
  int size=(2*radius+1)*(2*radius+1)-2*radius+1
  
  int x=1, y=1;    //initial position
  int dx=-1, dy=0; //iterative increment
  int c=2;         //number of positions per branch
  int d=0;         //directions: 0-left; 1-up; 2-right; 3-down
  
  int i=0
  index[i]=y*nx+x;

  while(i<size)
  {
    for(int j=0; j<c; j++)
    {
      x+=dx; y+=dy;
      if(x!=0) //do not include the current line in the index
      {
        index[i]=y*nx+x;
        i++;
      }
    }
    
    //start with the following branch
    if(d==2) 
    {
      x+=dx; y+=dy;
      index[i]=y*nx+x;
    }
    
    //change direction
    d=(d+1)%4;
    switch(d)
    {
      case 1: dx=-1; dy=0;  break;
      case 2: dx=0;  dy=-1; break;
      case 3: dx=1;  dy=0;  break;
      case 4: dx=0;  dy=1;  break;  
    }
    
    //increase the number of cells to traverse
    c++; 
  }
}


/**
  *
  * Function for non-maximum suppression
  *
**/
void non_maximum_suppression(
  float *I,            // input image
  int   radius,        // window radius
  int   nx,            // number of columns of the image
  int   ny,            // number of rows of the image
  int   verbose,       // activate verbose mode
  int   forensics      // activate forensics mode  
)
{
  int *mask     = new int[nx*ny];
  int *scanline = new int[nx];
  int *skip     = new int[nx];
  
  
  delete []mask;
  delete []scanline;
  delete []skip;
}


/**
  *
  * Function for computing Harris corners
  *
**/
void harris(
  float *I,            // input image
  std::vector<int> &x, // output selected points (x coordinates)
  std::vector<int> &y, // output selected points (y coordinates)
  float k,             // Harris constant for the ....function
  float sigma_i,       // standard deviation for smoothing (image denoising)    
  float sigma_n,       // standard deviation for smoothing (pixel neighbourhood)
  int   radius,        // radius of the autocorralation matrix
  float percentage,    
  int   nobel_measure,
  int   nx,            // number of columns of the image
  int   ny,            // number of rows of the image
  int   verbose,       // activate verbose mode
  int   forensics      // activate forensics mode  
)
{
  int size=nx*ny;
  float *Ix=new float[size];
  float *Iy=new float[size];
  float *Mc=new float[size];

  struct timeval start, end;
  
  //smooth the original image
  if (verbose)
  {
     printf("Harris corner detection:\n");
     printf(" 1.Convolving image with a Gaussian function\n");
    
     gettimeofday(&start, NULL);    
  }
  
  gaussian(I, nx, ny, sigma_i);
  
  
  if (verbose)
  {  
     gettimeofday(&end, NULL);
     printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
  }

  
  //compute the gradient of the image
  if (verbose) 
  {
     printf(" 2.Computing the gradient of the image\n");
  
     gettimeofday(&start, NULL);
  }
  
  gradient(I, Ix, Iy, nx, ny);
  
  
  if (verbose) 
  {  
     gettimeofday(&end, NULL);
     printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
  }

  if (forensics)
  {
     printf("  -Saving Is.png, Ix.png, Iy.png\n");
  
     char name1[200]="Is.png";
     char name2[200]="Ix.png";
     char name3[200]="Iy.png";
     iio_save_image_float_vec(name1, I, nx, ny, 1);
     iio_save_image_float_vec(name2, Ix, nx, ny, 1);
     iio_save_image_float_vec(name3, Iy, nx, ny, 1);
  }
  
  //compute the Harris function in each pixel
  if (verbose) 
  {
     printf(" 3.Computing the Harris' Mc function \n");

     gettimeofday(&start, NULL);
  }
  
  float *A = new float[size];
  float *B = new float[size];
  float *C = new float[size];

  #pragma omp parallel for
  for (int i=0;i<size;i++)
  {
     A[i] = Ix[i]*Ix[i];
     B[i] = Ix[i]*Iy[i];
     C[i] = Iy[i]*Iy[i];
  }
  
  gaussian(A, nx, ny, sigma_n);
  gaussian(B, nx, ny, sigma_n);
  gaussian(C, nx, ny, sigma_n);
  
  
  if (verbose) 
  {  
     gettimeofday(&end, NULL);
     printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
     
     gettimeofday(&start, NULL);     
  }


  float max = FLT_MIN;
  float min = FLT_MAX;
  
  if (nobel_measure == 0)
     harris_measure( A, B, C, Mc, nx, ny, k, max, min );

  if (nobel_measure == 1)
     shi_tomasi_measure( A, B, C, Mc, nx, ny, 0.0 );

  if (nobel_measure == 2)
     zseliski_measure( A, B, C, Mc, nx, ny, max, min );

  
  
  if (verbose) 
  {  
     gettimeofday(&end, NULL);
     printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
     
     printf("  -Mc max=%f, Mc min=%f\n",max, min);     
  }    


  // Non-maximum suppression
  if (verbose) 
  {
     printf(" 4.Non-maximum suppression\n");
     gettimeofday(&start, NULL);     
  }
  
  int *local_max = new int[nx*ny]();
  for (int i=radius+1;i<ny-radius-1;i++)
    for (int j=radius+1;j<nx-radius-1;j++)
        if(
	Mc[i*nx+j]>Mc[(i-1)*nx+j-1] && Mc[i*nx+j]>Mc[(i-1)*nx+j] &&
	Mc[i*nx+j]>Mc[(i-1)*nx+j+1] && Mc[i*nx+j]>Mc[i*nx+j-1]   &&
	Mc[i*nx+j]>Mc[i*nx+j+1] && Mc[i*nx+j]>Mc[(i+1)*nx+j-1]   &&
	Mc[i*nx+j]>Mc[(i+1)*nx+j] && Mc[i*nx+j]>Mc[(i+1)*nx+j+1] &&
	Mc[i*nx+j]>0
      )
	local_max[i*nx+j]=1;

  
  if (verbose) 
  {
     gettimeofday(&end, NULL);
     printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
  }

  //select the points depending on a percentage of the maximum and local maxima
  if (verbose)
  {
     printf(
        " 5.Selecting corner points by non-maximum suppression "
        "and %f of maximum\n", percentage
     );
  
     gettimeofday(&start, NULL);
  }

  x.reserve(1000);
  y.reserve(1000);
  for(int i=radius+1;i<ny-radius-1;i++)
    for(int j=radius+1;j<nx-radius-1;j++)
      if(local_max[i*nx+j] && Mc[i*nx+j]>max*percentage)
      {
        x.push_back(j);
        y.push_back(i);
      }

  if (verbose)
  {      
     gettimeofday(&end, NULL);
     printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
  }

      
  if(forensics)
  {
    printf("  -Number of detected corner points: %ld\n", x.size());
    printf("  -Saving Mc.png \n");
    for(int i=0;i<nx*ny;i++) Mc[i]=255*(Mc[i]-min)/(max-min); 
    char name[200]="Mc.png";
    iio_save_image_float_vec(name, Mc, nx, ny, 1);

    printf("  -Saving selected_Mc.png \n");
    for(int i=0;i<nx*ny;i++) if(local_max[i]==0) Mc[i]=0; 
    char name1[200]="selected_Mc.png";
    iio_save_image_float_vec(name1, Mc, nx, ny, 1);
  }
  
  delete []local_max;
  delete []Mc;
  delete []Ix;
  delete []Iy;
  delete []A;
  delete []B;
  delete []C;
}
