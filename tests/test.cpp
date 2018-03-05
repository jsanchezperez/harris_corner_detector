#include<cstdio>
#include<vector>
#include<cmath>
#include <random>

#include "../src/harris.h"
#include "../src/gaussian.h"
#include "../src/gradient.h"
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
   for(int i=0; i<points1.size(); i++)
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
   for(int i=0; i<points2.size(); i++)
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
      for(int i=0; i<s1.size(); i++)
      {
        int j=0;
        
        //search the point in the other set
        while(j<s2.size() && distance(s1[i], s2[j])>=TOL*TOL) 
          j++;

        //test if we found it
        if(j<s2.size()) N++;
      }
   }
   else
   {
      for(int i=0; i<s2.size(); i++)
      {
        int j=0;
        
        //search the point in the other set
        while(j<s1.size() && distance(s1[j], s2[i])>=TOL*TOL) 
          j++;
        
        //test if we found it
        if(j<s1.size()) N++;
      }
   }

   delete []p_1;

   float den = min(s1.size(), s2.size());

   //printf("We have found %d out of %lu (%d) \n", N, s2.size(), s1.size());

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
 *  Draw the Harris' corners
 *
 */
void draw_points(
  float *I, 
  std::vector<harris_corner> &corners,
  int strategy,
  int cells,
  int nx, 
  int ny, 
  int nz,
  int radius
)
{
  if(strategy==DISTRIBUTED_N_CORNERS)
  {
    //draw cells limits
    for(int i=0; i<cells; i++)
    {
      int cellx=cells, celly=cells;
      if(cellx>nx) cellx=nx;
      if(celly>ny) celly=ny;

      float Dx=(float)nx/cellx;
      float dx=Dx;
      while(dx<nx)
      {
        if(nz>=3)
          for(int y=0;y<ny;y++)
          {
            I[(y*nx+(int)dx)*nz]=0;
            I[(y*nx+(int)dx)*nz+1]=0;
            I[(y*nx+(int)dx)*nz+2]=0;
          }
        else
          for(int y=0;y<ny;y++)
            I[y*nx+(int)dx]=0;  
        dx+=Dx;
      }
    
      float Dy=(float)ny/celly;
      float dy=Dy;
      while(dy<ny)
      {
        if(nz>=3)
          for(int x=0;x<nx;x++)
          {
            I[((int)dy*nx+x)*nz]=0;
            I[((int)dy*nx+x)*nz+1]=0;
            I[((int)dy*nx+x)*nz+2]=0;
          }    
        else
          for(int x=0;x<nx;x++)
            I[(int)dy*nx+x]=0;
        dy+=Dy;
      }
    }
  }

  int dif=128;
  //draw a cross for each corner
  for(unsigned int i=0;i<corners.size();i++)
  {
    int x=corners[i].x;
    int y=corners[i].y;

    int x0=(x-radius<0)?0: x-radius;
    int x1=(x+radius>=nx)?nx-1: x+radius;
    int y0=(y-radius<0)?0: y-radius;
    int y1=(y+radius>=ny)?ny-1: y+radius;

    if(nz>=3)
    {
      int color=(int)(I[(y*nx+x)*nz+2]+dif)%256;
      
      for(int j=x0;j<=x1;j++)
      {
        I[(y*nx+j)*nz]=0;
        I[(y*nx+j)*nz+1]=0;
        I[(y*nx+j)*nz+2]=color;
      }
     
      for(int j=y0;j<=y1;j++)
      {
        I[(j*nx+x)*nz]=0;
        I[(j*nx+x)*nz+1]=0;
        I[(j*nx+x)*nz+2]=color;
      }
    }
    else
    {
      float c=0;
      int N=0;
      for(int j=x0;j<=x1;j++, N++)
        c+=I[y*nx+j];
      
      for(int j=y0;j<=y1;j++, N++)
        c+=I[j*nx+x];
      c/=(float) N;
      int color=(int)(128+c)%256;
      for(int j=x0;j<=x1;j++)
        I[y*nx+j]=color;
      
      for(int j=y0;j<=y1;j++)
        I[j*nx+x]=color;
    }
  }  
}


/**
 *
 *  Execute tests for repeatability
 *
 */
void create_graphic(
  char   *image, 
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
 int nx, ny, nz;
 
 float *Ic=iio_read_image_float_vec(image, &nx, &ny, &nz);
 int   gaussian=FAST_GAUSSIAN;
 int   gradient=CENTRAL_DIFFERENCES;
 int   measure=HARRIS_MEASURE;
 float k=0.06;
 float sigma_d=1.0;
 float sigma_i=SIGMA_I;
 float threshold=10;
 int   strategy=ALL_CORNERS;
 int   cells=1;
 int   Nselect=2000;
 int   subpixel_precision=1;
 int   verbose=1;
 
 if(opencv) threshold*=5500;
 
 if (Ic!=NULL)
 {
   std::vector<harris_corner> corners0;
   float *I=new float[nx*ny];
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
          //float x=nx/2.-0.5;
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
          //VER SI INTRODUZCO ROTATION
          nparams=3;
          T[0]=T[1]=T[2]=0;
          add_noise(I, Ii, alpha[i], nx*ny);
          for(int l=0; l<nx*ny; l++)
          { 
            if(Ii[l]<0) Ii[l]=0;
            else if(Ii[l]>255) Ii[l]=255;
          } 
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
      
      
      char name[200];
      draw_points(Ii, cornersi, strategy, cells, nx, ny, 1, 2*sigma_i);
      sprintf(name, "tmp/image%.4d.png",i);
      iio_save_image_float_vec(name, Ii, nx, ny, 1);
      
      
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
   printf("Cerramos fichero\n");

   delete []I;
   delete []Ii;
   free (Ic);
 }
}



/**
 *
 *  Main function.
 *  This program executes all the tests for a single image
 *
 */
int main(int argc, char *argv[])
{
   int N = 101;
   int NTOL = 10;
   double TOL[10] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
   double inc;
   double alpha[N];
   int opencv = 1; 


   //affine-graphic for a range of xi-repeatability
   printf("Generating affine graphics:\n");
   double affine[N];
   inc=1/(N-2.0);
   affine[0]=0; 
   for(int i=1; i<N; i++)
      affine[i]=affine[i-1]+inc;
      
   create_graphic(argv[1], affine, N, TOL, NTOL, AFFINE, "affine.txt", opencv);

//return 0;   

   //illumination-graphic for a range of xi-repeatability
   printf("Generating illumination graphics:\n");
   double illumination[N];
   inc=(6-0.1)/(N-2.0);
   illumination[0]=1; illumination[1]=0.1;
   for(int i=2; i<N; i++)
      illumination[i]=illumination[i-1]+inc;
    
   create_graphic(
     argv[1], illumination, N, TOL, NTOL, ILLUMINATION, "light.txt", opencv
   );

//return 0;

   //noise-graphic for a range of xi-repeatability
   printf("Generating noise graphics:\n");
   double noise[N];
   inc=30./N;
   noise[0]=0;
   for(int i=1; i<N; i++)
     noise[i]=noise[i-1]+inc;
    
   create_graphic(argv[1], noise, N, TOL, NTOL, NOISE, "noise.txt", opencv);

//return 0;   

   //scale-graphic for a range of xi-repeatability
   //the first value must be removed
   printf("Generating scale graphics:\n");
   double scales[N];
   inc=(4-0.2)/(N-2.0);
   scales[0]=1; scales[1]=0.2;
   for(int i=2; i<N; i++)
     scales[i]=scales[i-1]+inc;
    
   create_graphic(argv[1], scales, N, TOL, NTOL, SCALE, "scale.txt", opencv);

//return 0;   //comparing central differences and sobel operator 

   printf("Generating gradient graphics:\n");
   inc=PI/(N-1.0);
   alpha[0]=0;
   for(int i=1; i<N; i++)
     alpha[i]=alpha[i-1]+inc;
    
   create_graphic(
     argv[1], alpha, N, TOL, NTOL, GRADIENT, 
     "central_gradient.txt", opencv, CENTRAL_DIFFERENCES
   );
   create_graphic(
     argv[1], alpha, N, TOL, NTOL, GRADIENT, 
     "sobel_gradient.txt", opencv, SOBEL_OPERATOR
   );
   
   
//return 0;

   //comparing harris, harris improved and fast harris improved 
   printf("Generating Gaussian graphics:\n");
   inc=PI/(N-1.0);
   alpha[0]=0;
   for(int i=1; i<N; i++)
     alpha[i]=alpha[i-1]+inc;
    
  create_graphic(
     argv[1], alpha, N, TOL, NTOL, GAUSSIAN, 
     "std_gaussian.txt", opencv, STD_GAUSSIAN
   );
   create_graphic(
     argv[1], alpha, N, TOL, NTOL, GAUSSIAN, 
     "fast_gaussian.txt", opencv, FAST_GAUSSIAN
   );
   create_graphic(
     argv[1], alpha, N, TOL, NTOL, GAUSSIAN, 
     "no_gaussian.txt", opencv, NO_GAUSSIAN
   );
   
   
//return 0;

   //rotation-graphic for a range of xi-repeatability
   printf("Generating rotation graphics:\n");
   inc=PI/(N-1.0);
   alpha[0]=0;
   for(int i=1; i<N; i++)
     alpha[i]=alpha[i-1]+inc;
    
   create_graphic(
     argv[1], alpha, N, TOL, NTOL, ROTATION, "rotation.txt", opencv
   );

//return 0;

}


