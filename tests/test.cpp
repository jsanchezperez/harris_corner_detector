// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2018, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// All rights reserved.


#include<cstdio>
#include<vector>
#include<cmath>
#include <random>
#include <omp.h>

#include "../src/harris.h"
#include "../src/gaussian.h"
#include "../src/gradient.h"
#include "../src/interpolation.h"
#include "bicubic_interpolation.h"
#include "harris_opencv.h"

extern "C"
{
#include "../src/iio.h"
}

#define PI 3.141592653589793238462
#define ROTATION 0
#define SCALE 1  
#define AFFINE 2 
#define ILLUMINATION 3 
#define NOISE 4
#define GAUSSIAN 5
#define GRADIENT 6
#define SUBPIXEL 7

#define SIGMA_D 1
#define SIGMA_I 2.5
#define RADIUS 2*SIGMA_I

using namespace std;



/**
  *
  *  Function to add Gaussian noise to an image
  * 
**/
void add_noise(
  float *u,    //input image
  float *v,    //output noisy image
  float sigma, //Gaussian standard deviation
  int   size   //number of pixels
) 
{  
  std::default_random_engine generator;
  std::normal_distribution<float> dist(0, sigma);

  #pragma omp parallel for
  for (int i=0; i< size; i++) 
    v[i] =  u[i] + dist(generator);
}


float distance(const harris_corner &p1, const harris_corner &p2)
{
  
  float dx=p1.x-p2.x;
  float dy=p1.y-p2.y;
  
  return dx*dx+dy*dy;
}

/**
 *
 *  Function to compute the repeatability rate between two feature sets
 *  Schmid et al. (2000), Evaluation of Interest Point Detectors
 *
 */
float repetability_rate(
  vector<harris_corner> &points1, 
  vector<harris_corner> &points2, 
  float *p, 
  int   nparams,
  float TOL,
  int   nx,
  int   ny
)
{
   vector<harris_corner> s1; 
   vector<harris_corner> s2; 
   
   float *p_1 = new float[nparams];
   
   inverse_transform(p, p_1, nparams);

   //select points in common scene parts for first set
   for(unsigned int i=0; i<points1.size(); i++)
   {
     float xp, yp;
     project(points1[i].x, points1[i].y, p_1, xp, yp, nparams);
     
     if(xp>RADIUS && xp<nx-RADIUS && yp>RADIUS && yp<ny-RADIUS)
     {
       harris_corner p;
       p.x=xp; p.y=yp;
       s1.push_back(p);
     }
   }

   //select points in common scene parts for second set
   for(unsigned int i=0; i<points2.size(); i++)
   {
     float xp, yp;
     project(points2[i].x, points2[i].y, p, xp, yp, nparams);
     
     if(xp>RADIUS && xp<nx-RADIUS && yp>RADIUS && yp<ny-RADIUS)
       s2.push_back(points2[i]);
   }

   //compute xi-repeatability using the smaller set
   int N=0;
   if(s1.size()<s2.size())
   {
      for(unsigned int i=0; i<s1.size(); i++)
      {
        unsigned int j=0;
        
        //search the point in the other set
        while(j<s2.size() && distance(s1[i], s2[j])>=TOL*TOL) 
          j++;

        //check if we found it
        if(j<s2.size()) N++;
      }
   }
   else
   {
      for(unsigned int i=0; i<s2.size(); i++)
      {
        unsigned int j=0;
        
        //search the point in the other set
        while(j<s1.size() && distance(s1[j], s2[i])>=TOL*TOL) 
          j++;
        
        //check if we found it
        if(j<s1.size()) N++;
      }
   }

   delete []p_1;

   float den = min(s1.size(), s2.size());

   if(den==0) 
     return 0;
   else
     return (float) (N/den);
}


/**
  *
  *  Function for converting an rgb image to grayscale levels
  * 
**/
void rgb2gray(
  float *rgb,  //input color image
  float *gray, //output grayscale image
  int   nx,     //number of columns
  int   ny,     //number of rows
  int   nz      //number of channels
)
{
  #pragma omp parallel for
  for(int i=0;i<nx*ny;i++)
    gray[i]=(0.2989*rgb[i*nz]+0.5870*rgb[i*nz+1]+0.1140*rgb[i*nz+2]);
}



/**
 *
 *  Execute tests for repeatability
 *
 */
void create_graphic(
  float  *Ic,
  int    nx, 
  int    ny, 
  int    nz,
  double *alpha, 
  int    Nalpha, 
  double *TOL, 
  int    NTOL,
  int    type,
  const char *file_name,
  int    opencv=0,
  int    parameter=0
)
{ 
 int   gaussian=FAST_GAUSSIAN;
 int   gradient=CENTRAL_DIFFERENCES;
 int   measure=HARRIS_MEASURE;
 float k=0.06;
 float sigma_d=SIGMA_D;
 float sigma_i=SIGMA_I;
 float threshold=130;
 int   strategy=ALL_CORNERS;
 int   cells=1;
 int   Nselect=1500;
 int   subpixel_precision=1;
 int   verbose=1;
 
 if(opencv) threshold*=5500;
 
 if (Ic!=NULL)
 {
   std::vector<harris_corner> corners0;
   float *I =new float[nx*ny];
   float *Ii=new float[nx*ny];

   if(nz>1)
     rgb2gray(Ic, I, nx, ny, nz);
   else
     for(int i=0;i<nx*ny;i++)
       I[i]=Ic[i];
   
   FILE *file;
   file=fopen(file_name, "w");

   for(int i=0; i<Nalpha; i++)
   {
      std::vector<harris_corner> cornersi;
      int nparams=8;
      float T[8];

      switch(type)
      {
        case GAUSSIAN: {
          //rotate image with Euclidean transform
          gaussian=parameter;
          nparams=3;
          float x=nx/2.-0.5;
          float y=ny/2.-0.5;
          T[2]=alpha[i];
          T[0]=-cos(T[2])*x+sin(T[2])*y+x;
          T[1]=-sin(T[2])*x-cos(T[2])*y+y;
          bicubic_interpolation(I, Ii, T, nparams, nx, ny, 1);
          break;
        }  
        case GRADIENT: {
          //rotate image with Euclidean transform
          gradient=parameter;
          nparams=3;
          float x=nx/2.-0.5;
          float y=ny/2.-0.5;
          T[2]=alpha[i];
          T[0]=-cos(T[2])*x+sin(T[2])*y+x;
          T[1]=-sin(T[2])*x-cos(T[2])*y+y;
          bicubic_interpolation(I, Ii, T, nparams, nx, ny, 1);
          break;
        }  
        case ROTATION: {
          //rotate image with Euclidean transform
          nparams=3;
          float x=nx/2.-0.5;
          float y=ny/2.-0.5;
          T[2]=alpha[i];
          T[0]=-cos(T[2])*x+sin(T[2])*y+x;
          T[1]=-sin(T[2])*x-cos(T[2])*y+y;
          bicubic_interpolation(I, Ii, T, nparams, nx, ny, 1);
          break;
        }  
        case SCALE: {
          //scale image with Similarity
          nparams=4; //similarity transform
          float x=nx/2.-0.5;
          float y=ny/2.-0.5;
          T[0]=x*(1-alpha[i]);
          T[1]=y*(1-alpha[i]);
          T[2]=alpha[i]-1;
          T[3]=0;
          bicubic_interpolation(I, Ii, T, nparams, nx, ny, 1);
          break;
        }  
        case AFFINE: {
          //scale image with Similarity
          nparams=6; //similarity transform
          float y=ny/2.-0.5;
          T[1]=T[2]=T[4]=T[5]=0;
          T[0]=-alpha[i]*y;
          T[3]=alpha[i];
          bicubic_interpolation(I, Ii, T, nparams, nx, ny, 1);
          break;
        }  
        case ILLUMINATION:
          nparams=3;
          T[0]=T[1]=T[2]=0;
          for(int l=0; l<nx*ny; l++) Ii[l]=std::min(alpha[i]*I[l], 255.0);
          break;
          
        case NOISE:
          nparams=3;
          T[0]=T[1]=T[2]=0;
          add_noise(I, Ii, alpha[i], nx*ny);
          for(int l=0; l<nx*ny; l++)
          { 
            if(Ii[l]<0) Ii[l]=0;
            else if(Ii[l]>255) Ii[l]=255;
          } 
          break;
          
        case SUBPIXEL:
          nparams=3;
          float x=nx/2.-0.5;
          float y=ny/2.-0.5;
          T[2]=alpha[i];
          T[0]=-cos(T[2])*x+sin(T[2])*y+x;
          T[1]=-sin(T[2])*x-cos(T[2])*y+y;
          bicubic_interpolation(I, Ii, T, nparams, nx, ny, 1);
          subpixel_precision=parameter;
          break;
      }

      //compute the reference set
      if(i==0) 
      {
        float *II=new float[nx*ny];
        
        if(type!=NOISE)
          add_noise(Ii, II, 2, nx*ny);
        else 
          for(int l=0; l<nx*ny; l++) II[l]=Ii[l];

        if(opencv)
          opencv_harris(
            II, corners0, k, sigma_d, sigma_i, 
            threshold, strategy, cells, Nselect, subpixel_precision, 
            nx, ny, verbose
          ); 
        else
          harris(
            II, corners0, gaussian, gradient, measure, k, sigma_d, sigma_i, 
            threshold, strategy, cells, Nselect, subpixel_precision, 
            nx, ny, verbose
          );
        delete []II;
      }

      //add noise for realistic simulation
      if(type!=NOISE)
        add_noise(Ii, Ii, 2, nx*ny);

      if(opencv)
        opencv_harris(
          Ii, cornersi, k, sigma_d, sigma_i, 
          threshold, strategy, cells, Nselect, subpixel_precision, 
          nx, ny, verbose
        ); 
      else
        harris(
          Ii, cornersi, gaussian, gradient, measure, k, sigma_d, sigma_i, 
          threshold, strategy, cells, Nselect, subpixel_precision, 
          nx, ny, verbose
        );
      
      fprintf(file,"%lf ", alpha[i]);
      printf("%lf ", alpha[i]);
      
      for(int j=0; j<NTOL; j++)
      {
        float result=
          repetability_rate(corners0, cornersi, T, nparams, TOL[j], nx, ny);
        
        fprintf(file,"%f ", result);
        printf("%f ", result);
      }
      fprintf(file,"\n");
      printf("\n");
   } 
   
   fclose(file);

   delete []I;
   delete []Ii;
 }
}



/**
 *
 *  Main function
 *  This program executes all the tests for a single image
 *
 */
int main(int argc, char *argv[])
{
   if(argc >= 2)
   {
      int N = 101;
      int NTOL = 10;
      double TOL[10] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
      double inc;
      double alpha[N];
      int opencv = 0;
      int nx, ny, nz;
      float *Ic=iio_read_image_float_vec(argv[1], &nx, &ny, &nz);

      //comparing harris, harris improved and fast harris improved 
      printf("Generating Gaussian graphics:\n");
      inc=PI/(N-1.0);
      alpha[0]=0;
      for(int i=1; i<N; i++)
        alpha[i]=alpha[i-1]+inc;
       
      create_graphic(
        Ic, nx, ny, nz, alpha, N, TOL, NTOL, GAUSSIAN, 
        "no_gaussian.txt", opencv, NO_GAUSSIAN
      );
      
      create_graphic(
        Ic, nx, ny, nz, alpha, N, TOL, NTOL, GAUSSIAN, 
        "std_gaussian.txt", opencv, STD_GAUSSIAN
      );

      create_graphic(
        Ic, nx, ny, nz, alpha, N, TOL, NTOL, GAUSSIAN, 
        "fast_gaussian.txt", opencv, FAST_GAUSSIAN
      );

      //comparing central differences and sobel operator 
      printf("Generating gradient graphics:\n");
      inc=PI/(N-1.0);
      alpha[0]=0;
      for(int i=1; i<N; i++)
        alpha[i]=alpha[i-1]+inc;
       
      create_graphic(
        Ic, nx, ny, nz, alpha, N, TOL, NTOL, GRADIENT, 
        "central_gradient.txt", opencv, CENTRAL_DIFFERENCES
      );
      
      create_graphic(
        Ic, nx, ny, nz, alpha, N, TOL, NTOL, GRADIENT, 
        "sobel_gradient.txt", opencv, SOBEL_OPERATOR
      );   
     
      
      //rotation-graphic for a range of xi-repeatability
      printf("Generating rotation graphics:\n");
      inc=PI/(N-1.0);
      alpha[0]=0;
      for(int i=1; i<N; i++)
        alpha[i]=alpha[i-1]+inc;
       
      create_graphic(
        Ic, nx, ny, nz, alpha, N, TOL, NTOL, ROTATION, "rotation.txt", opencv
      );

      //noise-graphic for a range of xi-repeatability
      printf("Generating noise graphics:\n");
      double noise[N];
      inc=30./N;
      noise[0]=0;
      for(int i=1; i<N; i++)
        noise[i]=noise[i-1]+inc;
       
      create_graphic(
         Ic, nx, ny, nz, noise, N, TOL, NTOL, NOISE, "noise.txt", opencv
      );

      //illumination-graphic for a range of xi-repeatability
      printf("Generating illumination graphics:\n");
      double illumination[N];
      inc=(6-0.1)/(N-2.0);
      illumination[0]=1; illumination[1]=0.1;
      for(int i=2; i<N; i++)
         illumination[i]=illumination[i-1]+inc;
       
      create_graphic(
        Ic, nx, ny, nz, illumination, N, TOL, NTOL, 
        ILLUMINATION, "light.txt", opencv
      );

      //affine-graphic for a range of xi-repeatability
      printf("Generating affine graphics:\n");
      double affine[N];
      inc=1/(N-2.0);
      affine[0]=0; 
      for(int i=1; i<N; i++)
         affine[i]=affine[i-1]+inc;
         
      create_graphic(
         Ic, nx, ny, nz, affine, N, TOL, NTOL, AFFINE, "affine.txt", opencv
      );

      //scale-graphic for a range of xi-repeatability
      //the first value must be removed
      printf("Generating scale graphics:\n");
      double scales[N];
      inc=(4-0.2)/(N-2.0);
      scales[0]=1; scales[1]=0.2;
      for(int i=2; i<N; i++)
        scales[i]=scales[i-1]+inc;
       
      create_graphic(
         Ic, nx, ny, nz, scales, N, TOL, NTOL, SCALE, "scale.txt", opencv
      );

      //subpixel-graphic for a range of xi-repeatability
      printf("Generating subpixel graphics:\n");
      inc=PI/(N-1.0);
      alpha[0]=0;
      for(int i=1; i<N; i++)
        alpha[i]=alpha[i-1]+inc;
       
      create_graphic(
        Ic, nx, ny, nz, alpha, N, TOL, NTOL, SUBPIXEL, 
        "subpixel_quadratic.txt", opencv, QUADRATIC_APPROXIMATION
      );
      
      create_graphic(
        Ic, nx, ny, nz, alpha, N, TOL, NTOL, SUBPIXEL, 
        "subpixel_quartic.txt", opencv, QUARTIC_INTERPOLATION
      );

      free(Ic);
   }
   
   return 0;
}


