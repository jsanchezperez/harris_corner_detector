#include "harris.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <sys/time.h> 


extern "C"
{
#include "iio.h"
}


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
  int precision=5 //defines the size of the window
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



void harris_slow(
  float *I,
  std::vector<int> &x,
  std::vector<int> &y,
  float k,
  float sigma,
  int   radius,
  float percentage,
  int   nobel_measure,
  int   nx,
  int   ny,
  int   verbose
)
{
  int size=nx*ny;
  float *Ix=new float[size];
  float *Iy=new float[size];
  float *Mc=new float[size]();
  
  //smooth the original image
  if(verbose)
  {
    printf("Harris corner detection:\n");
    printf(" 1.Convolving image with a Gaussian function\n");
  }
  
  gaussian(I, nx, ny, 1.0);
  
  //compute the gradient of the image
  if(verbose)
    printf(" 2.Computing the gradient of the image\n");
  
  gradient(I, Ix, Iy, nx, ny);

  if(verbose)
  {
    printf("  -Saving Is.png, Ix.png, Iy.png\n");
  
    char name1[200]="Is.png";
    char name2[200]="Ix.png";
    char name3[200]="Iy.png";
    iio_save_image_float_vec(name1, I, nx, ny, 1);
    iio_save_image_float_vec(name2, Ix, nx, ny, 1);
    iio_save_image_float_vec(name3, Iy, nx, ny, 1);
  }
  
  float max=-1E15;
  float min=+1E15;

  //compute the Harris function in each pixel
  if(verbose)
    printf(" 3.Computing the Harris' Mc function \n");

  int ksize=2*radius+1;
  float *w=new float[ksize*ksize];
  for(int i=-radius; i<=radius; i++)
    for(int j=-radius; j<=radius; j++)
    {
      const float x2=i*i+j*j;
      const float sigma2=sigma*sigma;
      w[(i+radius)*ksize+j+radius]=exp(-0.5*x2/sigma2);
    }

  #pragma omp parallel for
  for(int i=radius+1;i<ny-radius-1;i++)
  {
    for(int j=radius+1;j<nx-radius-1;j++)
    {
      float A1=0.0;
      float A2=0.0;
      float A3=0.0;
      
      int ii=0;
      //compute the autocorrelation matrix
      for(int k=i-radius;k<=i+radius;k++)
      {
	int jj=0;
        for(int l=j-radius;l<=j+radius;l++)
        {
	  int p=ii*ksize+jj;
          A1+=w[p]*Ix[k*nx+l]*Ix[k*nx+l];
          A2+=w[p]*Ix[k*nx+l]*Iy[k*nx+l];
          A3+=w[p]*Iy[k*nx+l]*Iy[k*nx+l];
	  jj++;
        }
        ii++;
      }
    
      float detA=A1*A3-A2*A2;
      float traceA=A1+A3;
      
      if(nobel_measure)
	//Noble's corner measure (k is not needed)
	Mc[i*nx+j]=2*detA/(traceA+0.0001);
      else
	//Harris's corner measure
	Mc[i*nx+j]=detA-k*traceA*traceA;  
      
      if(Mc[i*nx+j]>max)
	max=Mc[i*nx+j];
      if(Mc[i*nx+j]<min)
	min=Mc[i*nx+j];
    }
  }

  if(verbose)
    printf("  -Mc max=%f, Mc min=%f\n",max, min);

  //non-maximum suppression
  if(verbose)
    printf(" 4.Non-maximum suppression\n");
  
  int *local_max=new int[nx*ny]();
  for(int i=radius+1;i<ny-radius-1;i++)
    for(int j=radius+1;j<nx-radius-1;j++)
      if(
	Mc[i*nx+j]>Mc[(i-1)*nx+j-1] && Mc[i*nx+j]>Mc[(i-1)*nx+j] &&
	Mc[i*nx+j]>Mc[(i-1)*nx+j+1] && Mc[i*nx+j]>Mc[i*nx+j-1]   &&
	Mc[i*nx+j]>Mc[i*nx+j+1] && Mc[i*nx+j]>Mc[(i+1)*nx+j-1]   &&
	Mc[i*nx+j]>Mc[(i+1)*nx+j] && Mc[i*nx+j]>Mc[(i+1)*nx+j+1] &&
	Mc[i*nx+j]>0
      )
	local_max[i*nx+j]=1;


  //select the points depending on a percentage of the maximum and local maxima
  if(verbose)
    printf(
      " 5.Selecting corner points by non-maximum suppression "
      "and %f of maximum\n", percentage
    );
  
  x.reserve(1000);
  y.reserve(1000);
  for(int i=radius+1;i<ny-radius-1;i++)
    for(int j=radius+1;j<nx-radius-1;j++)
      if(local_max[i*nx+j] && Mc[i*nx+j]>max*percentage)
      {
	x.push_back(j);
	y.push_back(i);
      }

      
  if(verbose)
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
  delete []w;
}




/**
  *
  * Function for computing Harris corners
  *
**/
void harris(
  float *I,            //input image
  std::vector<int> &x, //output selected points (x coordinates)
  std::vector<int> &y, //output selected points (y coordinates)
  float k,             //Harris constant for the ....function
  float sigma,         //standard deviation for smoothing 
  int   radius,        //radius of the autocorralation matrix
  float percentage,    
  int   nobel_measure,
  int   nx,            //number of columns of the image
  int   ny,            //number of rows of the image
  int   verbose        //activate verbose mode
)
{
  int size=nx*ny;
  float *Ix=new float[size];
  float *Iy=new float[size];
  float *Mc=new float[size]();

  struct timeval start, end;
  
  //smooth the original image
  //if(verbose)
  {
    printf("Harris corner detection:\n");
    printf(" 1.Convolving image with a Gaussian function\n");
  }
  
  gettimeofday(&start, NULL);

  gaussian(I, nx, ny, 1.0);
  
  gettimeofday(&end, NULL);
  printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6);

  
  //compute the gradient of the image
  //if(verbose)
    printf(" 2.Computing the gradient of the image\n");
  
  gettimeofday(&start, NULL);
  
  gradient(I, Ix, Iy, nx, ny);
  
  gettimeofday(&end, NULL);
  printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6);

  if(verbose)
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
//  if(verbose)
    printf(" 3.Computing the Harris' Mc function \n");

  gettimeofday(&start, NULL);
  
  float *A1=new float[size];
  float *A2=new float[size];
  float *A3=new float[size];

  #pragma omp parallel for
  for(int i=0;i<size;i++)
  {
    A1[i]=Ix[i]*Ix[i];
    A2[i]=Ix[i]*Iy[i];
    A3[i]=Iy[i]*Iy[i];
  }
  
  gaussian(A1, nx, ny, sigma);
  gaussian(A2, nx, ny, sigma);
  gaussian(A3, nx, ny, sigma);
  
  gettimeofday(&end, NULL);
  printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6);
  
  float max=FLT_MIN;
  float min=FLT_MAX;

  gettimeofday(&start, NULL);

  for(int i=0;i<size;i++)
  {
    float detA=A1[i]*A3[i]-A2[i]*A2[i];
    float traceA=A1[i]+A3[i];
    
    //if(nobel_measure)
      //Noble's corner measure (k is not needed)
      Mc[i]=2*detA/(traceA+0.0001);
    
      //else
      //Harris's corner measure
      //Mc[i]=detA-k*traceA*traceA;  
    
    if(Mc[i]>max)
      max=Mc[i];
    if(Mc[i]<min)
      min=Mc[i];
  }

  gettimeofday(&end, NULL);
  printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6);


  if(verbose)
    printf("  -Mc max=%f, Mc min=%f\n",max, min);

  //non-maximum suppression
//  if(verbose)
    printf(" 4.Non-maximum suppression\n");
  
  gettimeofday(&start, NULL);
  int *local_max=new int[nx*ny]();
  for(int i=radius+1;i<ny-radius-1;i++)
    for(int j=radius+1;j<nx-radius-1;j++)
      if(
	Mc[i*nx+j]>Mc[(i-1)*nx+j-1] && Mc[i*nx+j]>Mc[(i-1)*nx+j] &&
	Mc[i*nx+j]>Mc[(i-1)*nx+j+1] && Mc[i*nx+j]>Mc[i*nx+j-1]   &&
	Mc[i*nx+j]>Mc[i*nx+j+1] && Mc[i*nx+j]>Mc[(i+1)*nx+j-1]   &&
	Mc[i*nx+j]>Mc[(i+1)*nx+j] && Mc[i*nx+j]>Mc[(i+1)*nx+j+1] &&
	Mc[i*nx+j]>0
      )
	local_max[i*nx+j]=1;

  gettimeofday(&end, NULL);
  printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6);

  //select the points depending on a percentage of the maximum and local maxima
  //if(verbose)
    printf(
      " 5.Selecting corner points by non-maximum suppression "
      "and %f of maximum\n", percentage
    );
  
  gettimeofday(&start, NULL);

  x.reserve(1000);
  y.reserve(1000);
  for(int i=radius+1;i<ny-radius-1;i++)
    for(int j=radius+1;j<nx-radius-1;j++)
      if(local_max[i*nx+j] && Mc[i*nx+j]>max*percentage)
      {
        x.push_back(j);
        y.push_back(i);
      }

   gettimeofday(&end, NULL);
   printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
          end.tv_usec - start.tv_usec) / 1.e6);

      
  if(verbose)
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
  delete []A1;
  delete []A2;
  delete []A3;
}
