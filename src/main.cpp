#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h> 

#include "harris.h"

extern "C"
{
#include "iio.h"
}


#define PAR_DEFAULT_K 0.06
#define PAR_DEFAULT_SIGMA 2.5
#define PAR_DEFAULT_PERCENTAGE 0.1
#define PAR_DEFAULT_RADIUS 5
#define PAR_DEFAULT_NOBEL_MEASURE 0
#define PAR_DEFAULT_VERBOSE 0


using namespace std;


/**
 *
 *  Print a help message 
 *
 */
void print_help(char *name)
{
  printf("\n  Usage: %s image [OPTIONS] \n\n",
          name);
  printf("  Harris corner detector:\n");
  printf("  'image' is an input image to detect features on.\n");
  printf("  -----------------------------------------------\n");
  printf("  OPTIONS:\n"); 
  printf("  --------\n");
  printf("   -o name  output image with detected corners \n");
  printf("   -f name  write points to file\n");
  printf("   -k N     Harris' K parameter\n");
  printf("              default value %f\n", PAR_DEFAULT_K);
  printf("   -s N     Gauss standard deviation for weighted windows\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA);
  printf("   -w N     window radius size\n");
  printf("              default value %d\n", PAR_DEFAULT_RADIUS);
  printf("   -p N     percentage with respect to the maximum for selecting points\n");
  printf("              default value %f\n", PAR_DEFAULT_PERCENTAGE);
  printf("   -n       use Nobel's measure (does not need parameter K) \n");
  printf("   -v       switch on verbose mode \n\n\n");
}



/**
 *
 *  Read command line parameters 
 *
 */
int read_parameters(
  int   argc, 
  char  *argv[], 
  char  **image,
  char  **out_image,
  char  **out_file,
  float &k,
  float &sigma,
  int   &radius,
  float &percentage,
  int   &nobel_measure,
  int   &verbose
)
{
  if (argc < 2){
    print_help(argv[0]); 
    return 0;
  }
  else{
    int i=1;
    *image=argv[i++];

    //assign default values to the parameters
    k=PAR_DEFAULT_K;
    sigma=PAR_DEFAULT_SIGMA;
    radius=PAR_DEFAULT_RADIUS;
    verbose=PAR_DEFAULT_VERBOSE;
    percentage=PAR_DEFAULT_PERCENTAGE;
    nobel_measure=PAR_DEFAULT_NOBEL_MEASURE;
    
    //read each parameter from the command line
    while(i<argc)
    {
      if(strcmp(argv[i],"-o")==0)
        if(i<argc-1)
          *out_image=argv[++i];

      if(strcmp(argv[i],"-f")==0)
        if(i<argc-1)
          *out_file=argv[++i];
      
      if(strcmp(argv[i],"-a")==0)
        if(i<argc-1)
          k=atof(argv[++i]);

      if(strcmp(argv[i],"-s")==0)
        if(i<argc-1)
          sigma=atof(argv[++i]);
        
      if(strcmp(argv[i],"-w")==0)
        if(i<argc-1)
          radius=atoi(argv[++i]);

      if(strcmp(argv[i],"-p")==0)
        if(i<argc-1)
          percentage=atof(argv[++i]);

      if(strcmp(argv[i],"-n")==0)
        nobel_measure=1;

      if(strcmp(argv[i],"-v")==0)
        verbose=1;
      
      i++;
    }

    //check parameter values
    if(k<=0)
      k=PAR_DEFAULT_K;
    if(sigma<0)
      sigma=PAR_DEFAULT_SIGMA;
    if(radius<1)
      radius=PAR_DEFAULT_RADIUS;
    if(percentage<0 || percentage>1)
      percentage=PAR_DEFAULT_PERCENTAGE;
  }

  return 1;
}



void draw_points(
  float *I, 
  vector<int> &x, 
  vector<int> &y, 
  int nx, 
  int nz,
  int radius
)
{
  //draw a cross for each corner
  if(nz>=3)
    for(unsigned int i=0;i<x.size();i++)
    {
      for(int j=x[i]-radius;j<=x[i]+radius;j++)
      {
	I[(y[i]*nx+j)*nz]=0;
	I[(y[i]*nx+j)*nz+1]=0;
	I[(y[i]*nx+j)*nz+2]=255;
      }
      
      for(int j=y[i]-radius;j<=y[i]+radius;j++)
      {
	I[(j*nx+x[i])*nz]=0;
	I[(j*nx+x[i])*nz+1]=0;
	I[(j*nx+x[i])*nz+2]=255;
      }
    }
  else
    for(unsigned int i=0;i<x.size();i++)
    {
      for(int j=x[i]-radius;j<=x[i]+radius;j++)
	I[y[i]*nx+j]=255;
      for(int j=y[i]-radius;j<=y[i]+radius;j++)
	I[j*nx+x[i]]=255;
    }
}




/**
  *
  *  Function to convert an rgb image to grayscale levels
  * 
**/
void rgb2gray(
  float *rgb, //input color image
  float *gray, //output grayscale image
  int size     //number of pixels
)
{
  #pragma omp parallel for
  for(int i=0;i<size;i++)
    gray[i]=(0.2989*rgb[i*3]+0.5870*rgb[i*3+1]+0.1140*rgb[i*3+2]);
}



    
int main(int argc, char *argv[]) 
{
  //parameters of the method
  char  *image, *out_image=NULL, *out_file=NULL;
  float k, sigma, percentage;
  int   radius, nobel_measure, verbose;
    
  //read the parameters from the console
  int result=read_parameters(
        argc, argv, &image, &out_image, &out_file, 
	k, sigma, radius, percentage, nobel_measure, verbose
      );
  
  if(result)
  {
    int nx, ny, nz;   
    float *Ic=iio_read_image_float_vec(image, &nx, &ny, &nz);
    
    if(verbose)
      printf(
	"Parameters:\n"
        "  Input image: '%s', Output image: '%s', Output corner file: %s\n"
        "  K: %f, Sigma: %f, Window radius: %d, Percentage: %f\n"
	"  Nobel's measure: %d, nx: %d, ny: %d, nz: %d\n",
        image, out_image, out_file, k, sigma, radius, percentage,
	nobel_measure, nx, ny, nz
      );
    
    
    if (Ic!=NULL)
    {
      vector<int> x,y;
      float *I=new float[nx*ny];
      
      if(nz>1)
	rgb2gray(Ic, I, nx*ny);
      else
	for(int i=0;i<nx*ny;i++)
	  I[i]=Ic[i];

      struct timeval start, end;

      //if(verbose) 
        gettimeofday(&start, NULL);

      //compute Harris' corners
      harris(
	I, x, y, k, sigma, radius, percentage, 
	nobel_measure, nx, ny, verbose
      );
      
      //if(verbose) 
      {
        gettimeofday(&end, NULL);
        float delay=((end.tv_sec-start.tv_sec)* 1000000u + 
                     end.tv_usec - start.tv_usec) / 1.e6; 
        printf("\n Time: %fs\n", delay);
      }
  
      if(out_image!=NULL)
      {
	draw_points(Ic, x, y, nx, nz, radius);
	iio_save_image_float_vec(out_image, Ic, nx, ny, nz);
      }
	
      if(out_file!=NULL)
      {
	FILE *fd=fopen(out_file,"w");
	for(unsigned int i=0;i<x.size();i++)
	  fprintf(fd, "%d %d\n", x[i],y[i]);
	fclose(fd);
      }

      delete []I;
      free (Ic);
    }
    else 
    {
      printf("Cannot read image %s\n", image);
      exit(EXIT_FAILURE);
    }
  }
  exit(EXIT_SUCCESS);
}