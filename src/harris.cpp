#include <stdio.h>
#include <math.h>
#include <float.h>
#include <sys/time.h> 
#include <algorithm>

extern "C"
{
#include "iio.h"
}

#include "harris.h"
#include "gradient.h"
#include "gaussian.h"
#include "interpolation.h"


using namespace std;


/**
  *
  * Overload less function to compare two corners
  *
**/
bool operator<(
  const harris_corner &c1, 
  const harris_corner &c2
) {
  return c1.Mc > c2.Mc;
}


/**
  *
  * Function for computing the Autocorrelation matrix
  *
**/
void compute_autocorrelation_matrix(
  float *Ix,   //gradient of the image
  float *Iy,   //gradient of the image
  float *A,    //upper-left coefficient of the Autocorrelation matrix
  float *B,    //symmetric coefficient of the Autocorrelation matrix
  float *C,    //bottom-right coefficient of the Autocorrelation matrix
  float sigma, //standard deviation for smoothing (pixel neighbourhood)
  int   nx,    //number of columns of the image
  int   ny     //number of rows of the image
)
{
#pragma omp parallel for
  for (int i=0;i<nx*ny;i++)
  {
     A[i] = Ix[i]*Ix[i];
     B[i] = Ix[i]*Iy[i];
     C[i] = Iy[i]*Iy[i];
  }

  gaussian(A, nx, ny, sigma);
  gaussian(B, nx, ny, sigma);
  gaussian(C, nx, ny, sigma);
} 


/**
  *
  * Function for computing Harris' discriminant function
  *
**/
void compute_discriminant_function(
  float *A,      //upper-left coefficient of the Autocorrelation matrix
  float *B,      //symmetric coefficient of the Autocorrelation matrix
  float *C,      //bottom-right coefficient of the Autocorrelation matrix
  float *Mc,     //Harris measure
  int   measure, //measure strategy
  int   nx,      //number of columns of the image
  int   ny,      //number of rows of the image
  float k,       //Harris coefficient for the measure function
  float &max,    //output max value
  float &min     //output min value
)
{
  int size = nx*ny;

  //compute the discriminant function following one strategy
  switch(measure) 
  {
    default: case HARRIS_MEASURE:
      for (int i=0; i<size; i++)
      {
        float detA   = A[i]*C[i] - B[i]*B[i];
        float traceA = A[i] + C[i];

        Mc[i] = detA - k*traceA*traceA;

        if(Mc[i]>max) max=Mc[i];
        if(Mc[i]<min) min=Mc[i];
      }

    case HARMONIC_MEAN_MEASURE: 
      for (int i=0; i<size; i++)
      {
        float detA  =A[i]*C[i]-B[i]*B[i];
        float traceA=A[i]+C[i];

        Mc[i]=2*detA/(traceA+0.0001);

        if(Mc[i]>max) max=Mc[i];
        if(Mc[i]<min) min=Mc[i];
      }
      break;

    case SHI_TOMASI_MEASURE:
      for (int i=0; i<size; i++)
      {
        float D = sqrt(A[i]*A[i]-2*A[i]*C[i]+4*B[i]*B[i]+C[i]*C[i]);
        float lmin = 0.5*(A[i]+C[i])-0.5*D;

        if(lmin>k)
          Mc[i]=lmin;
        else
          Mc[i]=0;

        if(Mc[i]>max) max=Mc[i];
        if(Mc[i]<min) min=Mc[i];
      }
     break;
  }
}


/**
  *
  * Function for computing the spiral order indexes
  * for non-maximum_suppression
  *
**/
// void spiral_order(
//   int *index,
//   int radius,
//   int nx,
//   int bruteforce=0
// )
// {
//   int size=(2*radius+1)*(2*radius+1)-2*radius-1;
//   
//   int x=1, y=1;    //initial position
//   int dx=-1, dy=0; //iterative increment
//   int c=2;         //number of positions per branch
//   int d=0;         //directions: 0-left; 1-up; 2-right; 3-down
//   
//   //the first position is the botton-right corner
//   index[0]=y*nx+x;
//   
//   int i=1;
//   while(i<size)
//   {
//     //process following line
//     for(int j=0; j<c && i<size; j++)
//     {
//       //next position in the line
//       x+=dx; y+=dy;
//       
//       //do not include the central line in the index
//       if(y!=0 || bruteforce) 
//       {
//         index[i]=y*nx+x;
//         i++;
//       }
//     }    
//         
//     //change direction
//     d=(d+1)%4;
//     switch(d)
//     {
//       case 0: dx=-1; dy=0;  break;
//       case 1: dx=0;  dy=-1; break;
//       case 2: dx=1;  dy=0;  break;
//       case 3: dx=0;  dy=1;  break;  
//     }
//     
//     //increase traverse
//     if(d==0 || d==2) 
//       c++;
//   }
// }


//#define BRUTEFORCE
//#define SPIRAL_TEST

/**
  *
  * Function for non-maximum suppression
  *
**/
void non_maximum_suppression(
  float *I,             // input image
  vector<harris_corner> &corners,  // Harris' corners
  int   radius,         // window radius
  int   nx,             // number of columns of the image
  int   ny,             // number of rows of the image
  int   verbose         // activate verbose mode
)
{

// #ifdef BRUTEFORCE
//   
//   for(int i=radius; i<ny-radius; i++)
//   {
//     for(int j=radius; j<nx-radius; j++)
//     {
//       int k=i-radius; 
//       bool found=false;
//       while(!found && k<=i+radius)
//       {
//         int l=j-radius;
//         while(!found && l<=j+radius)
//         {
//           if(I[k*nx+l]>I[i*nx+j])
//             found=true;
//           l++;
//         }
//         k++;
//       }
//       
//       //if we found a local maximum add to the list
//       if(!found) 
//         corners.push_back({(float)j,(float)i,I[i*nx+j],0});
//     }
//   }
// 
// #else

  int *skip  = new int[nx*ny]();
  int size   = (2*radius+1)*(2*radius+1)-2*radius-1;
  //int *index = new int[size];

  //create the spiral order index
 // spiral_order(index, radius, nx);

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
// #ifdef SPIRAL_TEST
//             //spiral order test
//             int s=0;
//             while(s<size && I[i*nx+j+index[s]]<I[i*nx+j])
//             {
//               //if(i*nx+j+index[s]>i*nx+j)
//               skip[i*nx+j+index[s]]=1;
//               s++;	       
//             }
// 
//             //new local maximum found
//             if(s>=size) 
//               corners.push_back({(float)j,(float)i,I[i*nx+j],0});
// #else
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
              corners.push_back({(float)j,(float)i,I[i*nx+j],0});
// #endif
	  }
          j++;
        }
        else j=p1;
      }
    }
  }
  
  delete []skip;
  //delete []index;

//#endif
}


/**
  *
  * Function for computing subpixel precision of maxima
  *
**/
void subpixel_precision(
  float *Mc,  // discriminant function
  vector<harris_corner> &corners, // selected corners
  int nx
)
{
  for(unsigned int i=0; i<corners.size(); i++)
  {
    const float x=corners[i].x;
    const float y=corners[i].y;
   
    int mx=x-1;
    int dx=x+1;
    int my=y-1;
    int dy=y+1;
    
    float M[9];
    M[0]=Mc[my*nx+mx];
    M[1]=Mc[my*nx+(int)x];
    M[2]=Mc[my*nx+dx];
    M[3]=Mc[(int)y*nx+mx];
    M[4]=Mc[(int)y*nx+(int)x];
    M[5]=Mc[(int)y*nx+dx];
    M[6]=Mc[dy*nx+mx];
    M[7]=Mc[dy*nx+(int)x];
    M[8]=Mc[dy*nx+dx];
    
    corners[i].Mcint=corners[i].Mc;
    maximum_interpolation(M, corners[i]);
  }   
}


/**
  *
  * Function for selecting the output corners
  *
**/    
void select_output_corners(
  vector<harris_corner> &corners, // output selected corners
  int select_strategy, // strategy for the output corners
  int cells,           // number of regions in the image for distributed output
  int Nselect,         // number of output corners
  int nx,              // number of columns of the image
  int ny               // number of rows of the image
)
{
  switch(select_strategy) 
  {
    default: case ALL_CORNERS: 
      break;

    case ALL_CORNERS_SORTED:
      sort(corners.begin(), corners.end());
      break;

    case N_CORNERS:
      sort(corners.begin(), corners.end());
      corners.erase(corners.begin()+Nselect, corners.end());
      break;

    case DISTRIBUTED_N_CORNERS:
      int size=cells*cells;
      int Ncell=Nselect/size;
      vector<vector<harris_corner>> cell_corners(size);

      //distribute corners in the cells
      int Dx=nx/cells;
      int Dy=ny/cells;
      for(unsigned int i=0; i<corners.size(); i++)
      {
        int px=corners[i].x/Dx;
        int py=corners[i].y/Dy;

        cell_corners[py*cells+px].push_back(corners[i]);
      }

      //sort the corners in each cell
      for(int i=0; i<size; i++)
        sort(cell_corners[i].begin(), cell_corners[i].end());

      //copy the Ncell first corners to the output array
      corners.resize(0);
      for(int i=0; i<size; i++)
        corners.insert(
          corners.end(), cell_corners[i].begin(), 
          cell_corners[i].begin()+Ncell
        );
      break;
  }
}


/**
  *
  * Main function for computing Harris corners
  *
**/
void harris(
  float *I,         // input image
  vector<harris_corner> &corners, // output selected corners
  int   measure,    // measure for the discriminant function
  float k,          // Harris constant for the ....function
  float sigma_i,    // standard deviation for smoothing (image denoising)    
  float sigma_n,    // standard deviation for smoothing (pixel neighbourhood)
  int   radius,     // radius of the autocorralation matrix
  int   select_strategy, // strategy for the output corners
  int   cells,      // number of regions in the image for distributed output
  int   Nselect,    // number of output corners
  int   precision,  // enable subpixel precision
  int   nx,         // number of columns of the image
  int   ny,         // number of rows of the image
  int   verbose     // activate verbose mode
)
{
  int size=nx*ny;
  float *Ix=new float[size];
  float *Iy=new float[size];
  float *Mc=new float[size];

  struct timeval start, end;
  
  if (verbose)
  {
     printf("Harris corner detection:\n");
     printf(" 1.Convolving image with Gaussian (sigma=%f)\n", sigma_i);
    
     gettimeofday(&start, NULL);    
  }

  //smooth the original image to reduce noise
  gaussian(I, nx, ny, sigma_i);

  if (verbose)
  {  
     gettimeofday(&end, NULL);
     printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);

     printf(" 2.Computing the gradient of the image\n");
  
     gettimeofday(&start, NULL);
  }
  
  //compute the gradient of the image
  gradient(I, Ix, Iy, nx, ny);

  if (verbose) 
  {  
     gettimeofday(&end, NULL);
     printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);

     printf("  -Saving Is.png, Ix.png, Iy.png\n");

     char name1[200]="Is.png";
     char name2[200]="Ix.png";
     char name3[200]="Iy.png";
     iio_save_image_float_vec(name1, I, nx, ny, 1);
     iio_save_image_float_vec(name2, Ix, nx, ny, 1);
     iio_save_image_float_vec(name3, Iy, nx, ny, 1);

     printf(" 3.Computing the Harris' Mc function \n");

     gettimeofday(&start, NULL);
  }

  //variables for the Autocorrelation matrix
  float *A = new float[size];
  float *B = new float[size];
  float *C = new float[size];

  //compute the Autocorrelation matrix
  compute_autocorrelation_matrix(Ix, Iy, A, B, C, sigma_n, nx, ny);
  
  if (verbose) 
  {  
     gettimeofday(&end, NULL);
     printf("\n Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);

    printf(" 4.Computing one of the measures (0.Harris, 1.Shi-Tomasi, 2.Szeliski): %d\n", measure);
    gettimeofday(&start, NULL);     
  }

  float max  = FLT_MIN;
  float min  = FLT_MAX;

  //compute the discriminant function following one strategy 
  //Harris, Shi-Tomasi, Harmonic mean
  compute_discriminant_function(A, B, C, Mc, measure, nx, ny, k, max, min);

  if (verbose) 
  {
     gettimeofday(&end, NULL);
     printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
     
     printf("  -Mc max=%f, Mc min=%f\n",max, min);     
  }

  if (verbose) 
  {
     printf("\n 5.Non-maximum suppression\n");
     gettimeofday(&start, NULL);     
  }

  //select corners with non-maximum suppression 
  non_maximum_suppression(Mc, corners, radius, nx, ny, verbose);

  if (verbose) 
  {
     gettimeofday(&end, NULL);
     printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
  }

  if (verbose) 
  {
     printf("\n 6.Selecting output corners\n");
     gettimeofday(&start, NULL);     
  }

  //select output corners depending on the strategy
  //all corners; all corners sorted; N corners; distributed corners
  select_output_corners(corners, select_strategy, cells, Nselect, nx, ny);

  if (verbose) 
  {
     gettimeofday(&end, NULL);
     printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
  }

  if(precision)
  {
    if (verbose)
    {
       printf("\n 6.Computing subpixel precision\n");
       gettimeofday(&start, NULL);
    }

    //compute subpixel precision through quadratic interpolation
    subpixel_precision(Mc, corners, nx);

    if (verbose)
    {      
       gettimeofday(&end, NULL);
       printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);
    }
  }

  if(verbose)
  {
    printf("  -Number of detected corner points: %ld\n", corners.size());
    printf("  -Saving Mc.png \n");
    for(int i=0;i<nx*ny;i++) Mc[i]=255*(Mc[i]-min)/(max-min); 
    char name[200]="Mc.png";
    iio_save_image_float_vec(name, Mc, nx, ny, 1);

    printf("  -Saving selected_Mc.png \n");
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
