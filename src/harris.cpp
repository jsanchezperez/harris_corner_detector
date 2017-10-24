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
  printf("\n\n Shi-tomasi_measure: not implemented yet !!! \n\n");
}



/**
  *
  * Function for computing the spiral order indexes
  * for non-maximum_suppression
  *
**/
void spiral_order(
  int *index,
  int radius,
  int nx,
  int bruteforce=0
)
{
  int size=(2*radius+1)*(2*radius+1)-2*radius-1;
  
  int x=1, y=1;    //initial position
  int dx=-1, dy=0; //iterative increment
  int c=2;         //number of positions per branch
  int d=0;         //directions: 0-left; 1-up; 2-right; 3-down
  
  //the first position is the botton-right corner
  index[0]=y*nx+x;
  
  int i=1;
  while(i<size)
  {
    //process following line
    for(int j=0; j<c && i<size; j++)
    {
      //next position in the line
      x+=dx; y+=dy;
      
      //do not include the central line in the index
      if(y!=0 || bruteforce) 
      {
        index[i]=y*nx+x;
        i++;
      }
    }    
        
    //change direction
    d=(d+1)%4;
    switch(d)
    {
      case 0: dx=-1; dy=0;  break;
      case 1: dx=0;  dy=-1; break;
      case 2: dx=1;  dy=0;  break;
      case 3: dx=0;  dy=1;  break;  
    }
    
    //increase traverse
    if(d==0 || d==2) 
      c++;
  }
}


//#define BRUTEFORCE
//#define SPIRAL_TEST

/**
  *
  * Function for non-maximum suppression
  *
**/
void non_maximum_suppression(
  float *I,             // input image
  std::vector<float> &x,  // x position of maxima
  std::vector<float> &y,  // y position of maxima
  int   radius,         // window radius
  int   nx,             // number of columns of the image
  int   ny,             // number of rows of the image
  int   verbose         // activate verbose mode
)
{

#ifdef BRUTEFORCE
  
  for(int i=radius; i<ny-radius; i++)
  {
    for(int j=radius; j<nx-radius; j++)
    {
      int k=i-radius; 
      bool found=false;
      while(!found && k<=i+radius)
      {
        int l=j-radius;
        while(!found && l<=j+radius)
        {
          if(I[k*nx+l]>I[i*nx+j])
            found=true;
          l++;
        }
        k++;
      }
      
      //if we found a local maximum add to the list
      if(!found) 
      {
        x.push_back(j);
        y.push_back(i);
      }
    }
  }

#else

  int *skip  = new int[nx*ny]();
  int size   = (2*radius+1)*(2*radius+1)-2*radius-1;
  int *index = new int[size];

  //create the spiral order index
  spiral_order(index, radius, nx);

//#pragma omp parallel for
  for(int i=radius; i<ny-radius; i++)
  {
    int j=radius;
    
    //avoid the downhill at the beginning
    while(I[i*nx+j-1]>=I[i*nx+j] && j<nx-radius) j++;
      
    while(j<nx-radius)
    {
      //find the next peak 
      while((skip[i*nx+j] || I[i*nx+j+1]>=I[i*nx+j]) && j<nx-radius)
        j++;
      
      if(j<nx-radius)
      {
        int p1=j+2;

        //find a bigger value on the right
        while(I[i*nx+p1]<I[i*nx+j] && p1<=j+radius) 
        {
          skip[i*nx+p1]=1;
          p1++;
        }

        if(p1>j+radius)
        {  
          int p2=j-1;
  
          //find a bigger value on the left
          while(I[i*nx+p2]<=I[i*nx+j] && p2>=j-radius)
            p2--;
  
          if(p2<j-radius)
          {
#ifdef SPIRAL_TEST
            //spiral order test
            int s=0;
            while(s<size && I[i*nx+j+index[s]]<I[i*nx+j])
            {
              //if(i*nx+j+index[s]>i*nx+j)
              skip[i*nx+j+index[s]]=1;
              s++;	       
            }

            //new local maximum found
            if(s>=size) 
            {
              x.push_back(j);
              y.push_back(i);
            }
#else
            int k=i+radius; 
            bool found=false;
            while(!found && k>i)
            {
              int l=j+radius;
              while(!found && l>=j-radius)
              {
                if(I[k*nx+l]>I[i*nx+j])
                  found=true;
                else skip[k*nx+l]=1;
                l--;
              }
              k--;
            }
            
            k=i-radius; 
            while(!found && k<i)
            {
              int l=j-radius;
              while(!found && l<=j+radius)
              {
                if(I[k*nx+l]>=I[i*nx+j])
                  found=true;
                else skip[k*nx+l]=1;
                l++;
              }
              k++;
            }
            
            if(!found) 
            {
              x.push_back(j);
              y.push_back(i);
            }
#endif
	  }
          j++;
        }
        else j=p1;
      }
    }
  }
  
  delete []skip;
  delete []index;

#endif
}


/**
  *
  * Function for computing subpixel precision of maxima
  *
**/
void subpixel_precision(
  float *Mc,             // discriminant function
  std::vector<float> &x, // selected points (x coordinates)
  std::vector<float> &y  // selected points (y coordinates)
)
{
  /*for(int i=0; i<x.size(); i++)
  {
    float M[9];
    M[0]=
    //calculate the maximum of the quadratic interpolation
    maximum_interpolation(float *M, float x[i], float y[i])
  }*/
    
}

  
  
/**
  *
  * Function for computing Harris corners
  *
**/
void harris(
  float *I,            // input image
  std::vector<float> &x, // output selected points (x coordinates)
  std::vector<float> &y, // output selected points (y coordinates)
  float k,             // Harris constant for the ....function
  float sigma_i,       // standard deviation for smoothing (image denoising)    
  float sigma_n,       // standard deviation for smoothing (pixel neighbourhood)
  int   radius,        // radius of the autocorralation matrix
  float percentage,    
  int   nobel_measure,
  int   precision,
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
     
     
    printf(" 4.Computing one of the measures (0.Harris, 1.Shi-Tomasi, 2.Szeliski): %d\n", nobel_measure);
    gettimeofday(&start, NULL);     
  }

  float max = FLT_MIN;
  float min = FLT_MAX;
  
  //compute the discriminant function following one strategy
  if (nobel_measure == 0)
     harris_measure( A, B, C, Mc, nx, ny, k, max, min );

  if (nobel_measure == 1)
     shi_tomasi_measure( A, B, C, Mc, nx, ny, 0.0 );

  if (nobel_measure == 2)
     zseliski_measure( A, B, C, Mc, nx, ny, max, min );

  
  if (verbose) 
  {  
     gettimeofday(&end, NULL);
     printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
     
     printf("  -Mc max=%f, Mc min=%f\n",max, min);     
  }    


  // Non-maximum suppression
  if (verbose) 
  {
     printf("\n 5.Non-maximum suppression\n");
     gettimeofday(&start, NULL);     
  }
  
  x.reserve(1000);
  y.reserve(1000);
  non_maximum_suppression(Mc, x, y, radius, nx, ny, verbose);
  
  
  if (verbose) 
  {
     gettimeofday(&end, NULL);
     printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
  }


  //AQUI VA EL CODIGO PARA SELECCIONAR LOS PUNTOS (sort, ...)
  
 
  //compute subpixel precision through quadratic interpolation
  if(precision)
  {
    if (verbose)
    {
       printf("\n 6.Computing subpixel precision\n");
       gettimeofday(&start, NULL);
    }
    
    subpixel_precision(Mc, x, y);
    
    if (verbose)
    {      
       gettimeofday(&end, NULL);
       printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
    }    
    
  }
      
  if(forensics)
  {
    printf("  -Number of detected corner points: %ld\n", x.size());
    printf("  -Saving Mc.png \n");
    for(int i=0;i<nx*ny;i++) Mc[i]=255*(Mc[i]-min)/(max-min); 
    char name[200]="Mc.png";
    iio_save_image_float_vec(name, Mc, nx, ny, 1);

    printf("  -Saving selected_Mc.png \n");
    //for(int i=0;i<nx*ny;i++) if(local_max[i]==0) Mc[i]=0; 
    char name1[200]="selected_Mc.png";
    iio_save_image_float_vec(name1, Mc, nx, ny, 1);
  }
  
  delete []Mc;
  delete []Ix;
  delete []Iy;
  delete []A;
  delete []B;
  delete []C;
}
