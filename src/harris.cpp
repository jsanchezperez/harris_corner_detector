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
) 
{
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
  int   ny,    //number of rows of the image
  int   gauss  // type of Gaussian 

)
{
  #pragma omp parallel for
  for (int i=0;i<nx*ny;i++)
  {
     A[i] = Ix[i]*Ix[i];
     B[i] = Ix[i]*Iy[i];
     C[i] = Iy[i]*Iy[i];
  }
  
  if(gauss==NO_GAUSSIAN)
    gauss=FAST_GAUSSIAN;

  gaussian(A, nx, ny, sigma, gauss);
  gaussian(B, nx, ny, sigma, gauss);
  gaussian(C, nx, ny, sigma, gauss);
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
  float Th,      //threshold for eliminating low values
  float &max,    //output max value
  float &min     //output min value
)
{
  int size = nx*ny;

  //compute the discriminant function following one strategy
  switch(measure) 
  {
    default: case HARRIS_MEASURE:
      #pragma omp parallel for
      for (int i=0; i<size; i++)
      {
        float detA   = A[i]*C[i] - B[i]*B[i];
        float traceA = A[i] + C[i];

        Mc[i] = detA - k*traceA*traceA;

        if(Mc[i]>max) max=Mc[i];
        if(Mc[i]<min) min=Mc[i];
      }
      break;

    case HARMONIC_MEAN_MEASURE: 
      #pragma omp parallel for
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
      #pragma omp parallel for
      for (int i=0; i<size; i++)
      {
        float D = sqrt(A[i]*A[i]-2*A[i]*C[i]+4*B[i]*B[i]+C[i]*C[i]);
        float lmin = 0.5*(A[i]+C[i])-0.5*D;

        if(lmin>Th)
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
  * Function for non-maximum suppression
  *
**/
void non_maximum_suppression(
  float *D,             // input image
  vector<harris_corner> &corners, // Harris' corners
  int   radius,         // window radius
  int   nx,             // number of columns of the image
  int   ny              // number of rows of the image
)
{
  int *skip  = new int[nx*ny]();
  int size   = (2*radius+1)*(2*radius+1)-2*radius-1;
  
  //use an array for each row to allow parallel processing
  vector<vector<harris_corner> > corners_row(ny-2*radius);
 
  #pragma omp parallel for
  for(int i=radius; i<ny-radius; i++)
  {
    int j=radius;
    
    //avoid the downhill at the beginning
    while(D[i*nx+j-1]>=D[i*nx+j] && j<nx-radius) j++;
      
    while(j<nx-radius)
    {
      //find the next peak 
      while((skip[i*nx+j] || D[i*nx+j+1]>=D[i*nx+j]) && j<nx-radius)
        j++;
      
      if(j<nx-radius)
      {
        int p1=j+2;

        //find a bigger value on the right
        while(D[i*nx+p1]<D[i*nx+j] && p1<=j+radius) 
        {
          skip[i*nx+p1]=1;
          p1++;
        }

        //if not found
        if(p1>j+radius)
        {  
          int p2=j-1;
  
          //find a bigger value on the left
          while(D[i*nx+p2]<=D[i*nx+j] && p2>=j-radius)
            p2--;
  
          //if not found, test the 2D region
          if(p2<j-radius)
          {
            int k=i+radius; 
            bool found=false;

            //first test the bottom region (backwards)
            while(!found && k>i)
            {
              int l=j+radius;
              while(!found && l>=j-radius)
              {
                if(D[k*nx+l]>D[i*nx+j])
                  found=true;
                else skip[k*nx+l]=1;
                l--;
              }
              k--;
            }
            
            k=i-radius; 

            //then test the top region (forwards)
            while(!found && k<i)
            {
              int l=j-radius;
              while(!found && l<=j+radius)
              {
                if(D[k*nx+l]>=D[i*nx+j])
                  found=true;
                
                l++;
              }
              k++;
            }
            
            if(!found)
              //a new local maximum detected
              corners_row[i-radius].push_back({(float)j, (float)i, D[i*nx+j]});
          }
        }
        j=p1;
      }
    }
  }
  
  //copy row corners to the output list
  for(int i=0; i<ny-2*radius; i++)
   corners.insert(corners.end(), corners_row[i].begin(), corners_row[i].end());
  
  delete []skip;
}


/**
  *
  * Function for selecting the output corners
  *
**/    
void select_output_corners(
  vector<harris_corner> &corners, // output selected corners
  int strategy, // strategy for the output corners
  int cells,           // number of regions in the image for distributed output
  int Nselect,         // number of output corners
  int nx,              // number of columns of the image
  int ny               // number of rows of the image
)
{
  switch(strategy) 
  {
    default: case ALL_CORNERS: 
      break;

    case ALL_CORNERS_SORTED:
      sort(corners.begin(), corners.end());
      break;

    case N_CORNERS:
      sort(corners.begin(), corners.end());
      if(Nselect<(int)corners.size())
        corners.erase(corners.begin()+Nselect, corners.end());
      break;

    case DISTRIBUTED_N_CORNERS:
      int cellx=cells, celly=cells;
      if(cellx>nx) cellx=nx;
      if(celly>ny) celly=ny;

      int size=cellx*celly;
      int Ncell=Nselect/size;
      if(Ncell<1) Ncell=1;
      vector<vector<harris_corner>> cell_corners(size);
            
      float Dx=(float)nx/cellx;
      float Dy=(float)ny/celly;
      
      //distribute corners in the cells
      for(unsigned int i=0; i<corners.size(); i++)
      {
        int px=(float)corners[i].x/Dx;
        int py=(float)corners[i].y/Dy;
        cell_corners[(int)(py*cellx+px)].push_back(corners[i]);
      }

      //sort the corners in each cell
      for(int i=0; i<size; i++)
        sort(cell_corners[i].begin(), cell_corners[i].end());

      //copy the Ncell first corners to the output array
      corners.resize(0);
      for(int i=0; i<size; i++)
        if((int)cell_corners[i].size()>Ncell)
          corners.insert(
            corners.end(), cell_corners[i].begin(), 
            cell_corners[i].begin()+Ncell
          );
        else
          corners.insert(
            corners.end(), cell_corners[i].begin(), 
            cell_corners[i].end()
          );
      break;
  }
}


/**
  *
  * Function for computing subpixel precision of maxima
  *
**/
void compute_subpixel_precision(
  float *Mc, // discriminant function
  vector<harris_corner> &corners, // selected corners
  int nx,    // number of columns of the image
  int type   // type of interpolation (quadratic or quartic)
)
{
  #pragma omp parallel for
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
       
    //printf("c: %f, %f %f\n", corners[i].x, corners[i].y, corners[i].Mc);
    
    if(type==QUADRATIC_APPROXIMATION)
      quadratic_approximation(
        M, corners[i].x, corners[i].y, corners[i].Mc
      );
    else if(type==QUARTIC_INTERPOLATION)
      quartic_interpolation(
        M, corners[i].x, corners[i].y, corners[i].Mc
      );
    
    //printf("   %f, %f %f\n", corners[i].x, corners[i].y, corners[i].Mc);
    
  }   
}


/**
  *
  * Main function for computing Harris corners
  *
**/
void harris(
  float *I,        // input image
  vector<harris_corner> &corners, // output selected cornersç
  int   gauss,     // type of Gaussian 
  int   grad,      // type of gradient
  int   measure,   // measure for the discriminant function
  float k,         // Harris constant for the ....function
  float sigma_d,   // standard deviation for smoothing (image denoising)    
  float sigma_i,   // standard deviation for smoothing (pixel neighbourhood)
  float Th,        // threshold for eliminating low values
  int   strategy,  // strategy for the output corners
  int   cells,     // number of regions in the image for distributed output
  int   Nselect,   // number of output corners
  int   precision, // type of subpixel precision approximation
  int   nx,        // number of columns of the image
  int   ny,        // number of rows of the image
  int   verbose    // activate verbose mode
)
{
  int size=nx*ny;
  float *Ix=new float[size];
  float *Iy=new float[size];

  struct timeval start, end;
  
  if (verbose)
  {
    printf("\n\nHarris corner detection:\n");
    printf(" 1.Convolving image with Gaussian: \t");
    
    gettimeofday(&start, NULL);    
  }

  //smooth the original image to reduce noise
  gaussian(I, nx, ny, sigma_d, gauss);

  if (verbose)
  {  
    gettimeofday(&end, NULL);
    printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
           end.tv_usec - start.tv_usec) / 1.e6);

    printf(" 2.Computing the gradient: \t \t");
    gettimeofday(&start, NULL);    

  }
  
  //compute the gradient of the image
  gradient(I, Ix, Iy, nx, ny, grad);

  if (verbose) 
  {  
    gettimeofday(&end, NULL);
    printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);

//      printf("  -Saving Is.png, Ix.png, Iy.png\n");
// 
//      char name1[200]="Is.png";
//      char name2[200]="Ix.png";
//      char name3[200]="Iy.png";
//      iio_save_image_float_vec(name1, I, nx, ny, 1);
//      iio_save_image_float_vec(name2, Ix, nx, ny, 1);
//      iio_save_image_float_vec(name3, Iy, nx, ny, 1);

    printf(" 3.Computing the Autocorrelation: \t");

    gettimeofday(&start, NULL);
  }

  //variables for the Autocorrelation matrix
  float *A = new float[size];
  float *B = new float[size];
  float *C = new float[size];

  //compute the Autocorrelation matrix
  compute_autocorrelation_matrix(Ix, Iy, A, B, C, sigma_i, nx, ny, gauss);

  if (verbose) 
  {  
    gettimeofday(&end, NULL);
    printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
            end.tv_usec - start.tv_usec) / 1.e6);

    printf(" 4.Computing Harris function: \t\t");
    gettimeofday(&start, NULL);     
  }

  float max = FLT_MIN;
  float min = FLT_MAX;
  float *Mc = new float[size];

  //calculate the discriminant function (Harris, Shi-Tomasi, Harmonic mean)
  compute_discriminant_function(A, B, C, Mc, measure, nx, ny, k, Th, max, min);

  delete []Ix;
  delete []Iy;
  delete []A;
  delete []B;
  delete []C;

  if (verbose) 
  {
    gettimeofday(&end, NULL);
    printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
           end.tv_usec - start.tv_usec) / 1.e6);
     
    // printf("  -Mc max=%f, Mc min=%f\n",max, min);     
     
    printf(" 5.Apply threshold: \t\t \t");
    gettimeofday(&start, NULL);     
  }
  
  //threshold the discriminant function
  if(measure!=SHI_TOMASI_MEASURE)
    #pragma omp parallel for
    for(int i=0; i<size; i++)
      if(Mc[i]<Th) Mc[i]=0;
  
    
  //char name[200]= "harris_opencv.png";
  //iio_save_image_float_vec(name, Mc, nx, ny, 1);
      

      
  if (verbose) 
  {
    gettimeofday(&end, NULL);
    printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
           end.tv_usec - start.tv_usec) / 1.e6);
     
    printf(" 6.Non-maximum suppression:  \t\t");
    gettimeofday(&start, NULL);     
  }

  //apply non-maximum suppression to select salient corners
  int radius=2*sigma_i;
  non_maximum_suppression(Mc, corners, radius, nx, ny);

  if (verbose) 
  {
    gettimeofday(&end, NULL);
    printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
           end.tv_usec - start.tv_usec) / 1.e6);
  }

  if (verbose) 
  {
    printf(" 7.Selecting output corners:  \t\t");
    gettimeofday(&start, NULL);     
  }

  //select output corners depending on the strategy
  //all corners; all corners sorted; N corners; distributed corners
  select_output_corners(corners, strategy, cells, Nselect, nx, ny);

  if (verbose) 
  {
    gettimeofday(&end, NULL);
    printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
           end.tv_usec - start.tv_usec) / 1.e6);
  }

  if(precision==QUADRATIC_APPROXIMATION || precision==QUARTIC_INTERPOLATION)
  {
    if (verbose)
    {
      printf(" 8.Computing subpixel accuracy: \t");
      gettimeofday(&start, NULL);
    }

    //calculate subpixel precision
    compute_subpixel_precision(Mc, corners, nx, precision);

    if (verbose)
    {      
      gettimeofday(&end, NULL);
      printf("Time: %fs\n", ((end.tv_sec-start.tv_sec)* 1000000u + 
           end.tv_usec - start.tv_usec) / 1.e6);
    }
  }

  if(verbose)
  {
    printf(" * Number of corners detected: %ld\n", corners.size());
//     printf("  -Saving Mc.png \n");
//     for(int i=0;i<nx*ny;i++) Mc[i]=255*(Mc[i]-min)/(max-min); 
//     char name[200]="Mc.png";
//     iio_save_image_float_vec(name, Mc, nx, ny, 1);
// 
//     printf("  -Saving selected_Mc.png \n");
//     char name1[200]="selected_Mc.png";
//     iio_save_image_float_vec(name1, Mc, nx, ny, 1);
  }
  
  delete []Mc;
}
