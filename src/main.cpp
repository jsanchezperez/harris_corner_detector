#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h> 
#include <math.h>
#include <float.h>

#include "harris.h"

extern "C"
{
#include "iio.h"
}


#define PAR_DEFAULT_K 0.06
#define PAR_DEFAULT_SIGMA_I 1.0
#define PAR_DEFAULT_SIGMA_N 2.5
#define PAR_DEFAULT_RADIUS 5
#define PAR_DEFAULT_MEASURE HARRIS_MEASURE
#define PAR_DEFAULT_SELECT_STRATEGY N_CORNERS
#define PAR_DEFAULT_NSELECT 2000
#define PAR_DEFAULT_SUBPIXEL_PRECISION 0
#define PAR_DEFAULT_VERBOSE 0
#define PAR_DEFAULT_FORENSIC 0



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
  printf("   -m       choose measure: 0.Harris; 1.Shi-Tomasi; 2.Harmonic Mean\n"); 
  printf("   -k N     Harris' K parameter\n");
  printf("              default value %f\n", PAR_DEFAULT_K);
  printf("   -i N     Gauss standard deviation for image denoising\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_I);    
  printf("   -s N     Gauss standard deviation for weighted windows\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_N);
  printf("   -w N     window radius size\n");
  printf("              default value %d\n", PAR_DEFAULT_RADIUS);
  printf("   -q N     strategy for selecting the output corners\n");
  printf("              default value %d\n", PAR_DEFAULT_SELECT_STRATEGY);
  printf("   -N N     number of output corners\n");
  printf("              default value %d\n", PAR_DEFAULT_NSELECT);
  printf("   -p       subpixel precision through quadratic interpolation\n"); 
  printf("   -v       switch on verbose mode \n");
  printf("   -b       switch on forensic mode \n\n\n");  
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
  int   &measure,
  float &k,
  float &sigma_i,  
  float &sigma_n,
  int   &radius,
  int   &select_strategy,
  int   &Nselect,
  int   &subpixel_precision,  
  int   &verbose,
  int   &forensic 
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
    sigma_i=PAR_DEFAULT_SIGMA_I;    
    sigma_n=PAR_DEFAULT_SIGMA_N;
    measure=PAR_DEFAULT_MEASURE;
    radius=PAR_DEFAULT_RADIUS;
    select_strategy=PAR_DEFAULT_SELECT_STRATEGY;
    Nselect=PAR_DEFAULT_NSELECT;
    subpixel_precision=PAR_DEFAULT_SUBPIXEL_PRECISION;
    verbose=PAR_DEFAULT_VERBOSE;
    forensic=PAR_DEFAULT_FORENSIC;
    
    //read each parameter from the command line
    while(i<argc)
    {
      if(strcmp(argv[i],"-o")==0)
        if(i<argc-1)
          *out_image=argv[++i];

      if(strcmp(argv[i],"-f")==0)
        if(i<argc-1)
          *out_file=argv[++i];
      
      if(strcmp(argv[i],"-m")==0)
        measure=atoi(argv[++i]);

      if(strcmp(argv[i],"-k")==0)
        if(i<argc-1)
          k=atof(argv[++i]);

      if(strcmp(argv[i],"-i")==0)
        if(i<argc-1)
          sigma_i=atof(argv[++i]);        
        
      if(strcmp(argv[i],"-s")==0)
        if(i<argc-1)
          sigma_n=atof(argv[++i]);
        
      if(strcmp(argv[i],"-w")==0)
        if(i<argc-1)
          radius=atoi(argv[++i]);

      if(strcmp(argv[i],"-q")==0)
        if(i<argc-1)
          select_strategy=atoi(argv[++i]);
	
      if(strcmp(argv[i],"-N")==0)
        if(i<argc-1)
          Nselect=atoi(argv[++i]);

      if(strcmp(argv[i],"-p")==0)
        subpixel_precision=1;

      if(strcmp(argv[i],"-v")==0)
        verbose=1;
      
      if(strcmp(argv[i],"-b")==0)
        forensic=1;
      
      i++;
    }

    //check parameter values
    if (k<=0)        k       = PAR_DEFAULT_K;
    if (sigma_i<0)   sigma_i = PAR_DEFAULT_SIGMA_I;
    if (sigma_n<0)   sigma_n = PAR_DEFAULT_SIGMA_N;
    if (radius<1)    radius  = PAR_DEFAULT_RADIUS;    
    if (Nselect<1)   Nselect = PAR_DEFAULT_NSELECT;
  }

  return 1;
}



void draw_points(
  float *I, 
  std::vector<harris_corner> &corners, 
  int nx, 
  int nz,
  int radius
)
{
  float max=FLT_MIN;
  
  //find maximum of Harris' measure
  for(unsigned int i=0;i<corners.size();i++)
    if(corners[i].Mc>max) max=corners[i].Mc;
  
  //draw a cross for each corner
  if(nz>=3)
    for(unsigned int i=0;i<corners.size();i++)
    {
      const float x=corners[i].x;
      const float y=corners[i].y;
      const float C=std::max(50., 255.-255*pow(1-corners[i].Mc/max,5));
      
      for(int j=x-radius;j<=x+radius;j++)
      {
	I[((int)y*nx+j)*nz]=0;
	I[((int)y*nx+j)*nz+1]=0;
	I[((int)y*nx+j)*nz+2]=C;
      }
      
      for(int j=y-radius;j<=y+radius;j++)
      {
	I[(j*nx+(int)x)*nz]=0;
	I[(j*nx+(int)x)*nz+1]=0;
	I[(j*nx+(int)x)*nz+2]=C;
      }
    }
  else
    for(unsigned int i=0;i<corners.size();i++)
    {
      const float x=corners[i].x;
      const float y=corners[i].y;
      const float C=pow(0.3,1-corners[i].Mc/max);
      
      for(int j=x-radius;j<=x+radius;j++)
	I[(int)y*nx+j]=255*C;
      for(int j=y-radius;j<=y+radius;j++)
	I[j*nx+(int)x]=255*C;
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
  float k, sigma_i, sigma_n;
  int   select_strategy, Nselect, radius, measure, subpixel_precision;
  int   verbose, forensic;
    
  //read the parameters from the console
  int result=read_parameters(
        argc, argv, &image, &out_image, &out_file, 
        measure, k, sigma_i, sigma_n, radius, select_strategy, 
        Nselect, subpixel_precision, verbose, forensic
      );
  
  if(result)
  {
    int nx, ny, nz;   
    float *Ic=iio_read_image_float_vec(image, &nx, &ny, &nz);
    
    if(verbose)
      printf(
        "Parameters:\n"
        "  Input image: '%s', Output image: '%s', Output corner file: %s\n"
        "  measure: %d, K: %f, Sigma_i: %f, Sigma_n: %f, Window radius: %d, \n"
        "  Select strategy: %d, Nselect: %d, nx: %d, ny: %d, nz: %d\n",
        image, out_image, out_file,  measure, k, sigma_i, sigma_n, radius, 
	select_strategy, Nselect, nx, ny, nz
      );
    
    
    if (Ic!=NULL)
    {
      std::vector<harris_corner> corners;
      float *I=new float[nx*ny];
      
      if(nz>1)
	rgb2gray(Ic, I, nx*ny);
      else
	for(int i=0;i<nx*ny;i++)
	  I[i]=Ic[i];

      struct timeval start, end;

      if (verbose) 
        gettimeofday(&start, NULL);

      //compute Harris' corners
      harris(
        I, corners, measure, k, sigma_i, sigma_n, radius, select_strategy,
	Nselect, subpixel_precision, nx, ny, verbose, forensic
      );
      
      if (verbose) 
      {
        gettimeofday(&end, NULL);
        float delay=((end.tv_sec-start.tv_sec)* 1000000u + 
                     end.tv_usec - start.tv_usec) / 1.e6; 
        printf("\n Time: %fs\n", delay);
      }
  
      if(out_image!=NULL)
      {
	draw_points(Ic, corners, nx, nz, radius);
	iio_save_image_float_vec(out_image, Ic, nx, ny, nz);
      }
	
      if(out_file!=NULL)
      {
	FILE *fd=fopen(out_file,"w");
	for(unsigned int i=0;i<corners.size();i++)
	  fprintf(fd, "%f %f %f\n", corners[i].x, corners[i].y, corners[i].Mc);
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
