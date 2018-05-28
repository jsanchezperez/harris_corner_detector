#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h> 
#include <math.h>
#include <float.h>
#include <omp.h>

#include "harris.h"
#include "gaussian.h"
#include "gradient.h"
#include "interpolation.h"

extern "C"
{
#include "iio.h"
}

#define PAR_DEFAULT_K 0.06
#define PAR_DEFAULT_SIGMA_D 1.0
#define PAR_DEFAULT_SIGMA_I 2.5
#define PAR_DEFAULT_THRESHOLD 130
#define PAR_DEFAULT_GAUSSIAN FAST_GAUSSIAN
#define PAR_DEFAULT_GRADIENT CENTRAL_DIFFERENCES
#define PAR_DEFAULT_MEASURE HARRIS_MEASURE
#define PAR_DEFAULT_SELECT_STRATEGY ALL_CORNERS
#define PAR_DEFAULT_CELLS 3
#define PAR_DEFAULT_NSELECT 2000
#define PAR_DEFAULT_PRECISION QUADRATIC_APPROXIMATION
#define PAR_DEFAULT_VERBOSE 0

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
  printf("   -s N     choose smoothing: \n"); 
  printf("              0.precise Gaussian; 1.fast Gaussian; 2.no Gaussian\n"); 
  printf("              default value %d\n", PAR_DEFAULT_GAUSSIAN);
  printf("   -g N     choose gradient: \n"); 
  printf("              0.central differences; 1.Sobel operator\n"); 
  printf("              default value %d\n", PAR_DEFAULT_GRADIENT);
  printf("   -m N     choose measure: \n"); 
  printf("              0.Harris; 1.Shi-Tomasi; 2.Harmonic Mean\n"); 
  printf("              default value %d\n", PAR_DEFAULT_MEASURE);
  printf("   -k N     Harris' K parameter\n");
  printf("              default value %f\n", PAR_DEFAULT_K);
  printf("   -d N     Gaussian standard deviation for derivation\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_D);    
  printf("   -i N     Gaussian standard deviation for integration\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_I);
  printf("   -t N     threshold for eliminating low values\n");
  printf("              default value %d\n", PAR_DEFAULT_THRESHOLD);
  printf("   -q N     strategy for selecting the output corners:\n");
  printf("              0.all corners; 1.sort all corners;\n");
  printf("              2.N corners; 3.distributed N corners\n");
  printf("              default value %d\n", PAR_DEFAULT_SELECT_STRATEGY);
  printf("   -c N     regions for output corners (1x1, 2x2,...NxN):\n");
  printf("              default value %d\n", PAR_DEFAULT_CELLS);
  printf("   -n N     number of output corners\n");
  printf("              default value %d\n", PAR_DEFAULT_NSELECT);
  printf("   -p N     subpixel accuracy\n");
  printf("              0.no subpixel; 1.quadratic approximation;"
                      " 2.quartic interpolation\n");
  printf("              default value %d\n", PAR_DEFAULT_PRECISION);
  printf("   -v       switch on verbose mode \n");
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
  int   &gaussian,
  int   &gradient,
  int   &measure,
  float &k,
  float &sigma_d,  
  float &sigma_i,
  float &threshold,
  int   &strategy,
  int   &cells,
  int   &Nselect,
  int   &precision,  
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
    sigma_d=PAR_DEFAULT_SIGMA_D;    
    sigma_i=PAR_DEFAULT_SIGMA_I;
    gaussian=PAR_DEFAULT_GAUSSIAN;
    gradient=PAR_DEFAULT_GRADIENT;
    measure=PAR_DEFAULT_MEASURE;
    threshold=PAR_DEFAULT_THRESHOLD;
    strategy=PAR_DEFAULT_SELECT_STRATEGY;
    cells=PAR_DEFAULT_CELLS;
    Nselect=PAR_DEFAULT_NSELECT;
    precision=PAR_DEFAULT_PRECISION;
    verbose=PAR_DEFAULT_VERBOSE;
    
    //read each parameter from the command line
    while(i<argc)
    {
      if(strcmp(argv[i],"-o")==0)
        if(i<argc-1)
          *out_image=argv[++i];

      if(strcmp(argv[i],"-f")==0)
        if(i<argc-1)
          *out_file=argv[++i];
      
      if(strcmp(argv[i],"-s")==0)
        gaussian=atoi(argv[++i]);

      if(strcmp(argv[i],"-g")==0)
        gradient=atoi(argv[++i]);

      if(strcmp(argv[i],"-m")==0)
        measure=atoi(argv[++i]);

      if(strcmp(argv[i],"-k")==0)
        if(i<argc-1)
          k=atof(argv[++i]);

      if(strcmp(argv[i],"-d")==0)
        if(i<argc-1)
          sigma_d=atof(argv[++i]);        
        
      if(strcmp(argv[i],"-i")==0)
        if(i<argc-1)
          sigma_i=atof(argv[++i]);
        
      if(strcmp(argv[i],"-t")==0)
        if(i<argc-1)
          threshold=atof(argv[++i]);

      if(strcmp(argv[i],"-q")==0)
        if(i<argc-1)
          strategy=atoi(argv[++i]);

      if(strcmp(argv[i],"-c")==0)
        if(i<argc-1)
          cells=atoi(argv[++i]);

      if(strcmp(argv[i],"-n")==0)
        if(i<argc-1)
          Nselect=atoi(argv[++i]);

      if(strcmp(argv[i],"-p")==0)
        if(i<argc-1)
          precision=atoi(argv[++i]);

      if(strcmp(argv[i],"-v")==0)
        verbose=1;

      i++;
    }

    //check parameter values
    if(k<=0)      k       = PAR_DEFAULT_K;
    if(sigma_d<0) sigma_d = PAR_DEFAULT_SIGMA_D;
    if(sigma_i<0) sigma_i = PAR_DEFAULT_SIGMA_I;
    if(cells<1)   cells   = PAR_DEFAULT_CELLS;
    if(Nselect<1) Nselect = PAR_DEFAULT_NSELECT;
  }

  return 1;
}



/**
 *
 *  Draw the Harris' corners and division cells on the image
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
    #pragma omp parallel for
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

  //draw a cross for each corner
  //#pragma omp parallel for
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
      for(int j=x0;j<=x1;j++)
      {
        I[(y*nx+j)*nz]=0;
        I[(y*nx+j)*nz+1]=0;
        I[(y*nx+j)*nz+2]=255;
      }
    }
    else
      for(int j=x0;j<=x1;j++)
          I[y*nx+j]=255;
    
    if(nz>=3)
    {
      for(int j=y0;j<=y1;j++)
      {
        I[(j*nx+x)*nz]=0;
        I[(j*nx+x)*nz+1]=0;
        I[(j*nx+x)*nz+2]=255;
      }
    }
    else
        for(int j=y0;j<=y1;j++)
          I[j*nx+x]=255;
    }
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
 *  Main function
 *
 */
int main(int argc, char *argv[]) 
{      
  //parameters of the method
  char  *image, *out_image=NULL, *out_file=NULL;
  float k, sigma_d, sigma_i, threshold;
  int   gaussian, gradient, strategy, Nselect, measure;
  int   precision, cells, verbose;

  printf("Read parameters\n");
  //read the parameters from the console
  int result=read_parameters(
        argc, argv, &image, &out_image, &out_file, 
        gaussian, gradient, measure, k, sigma_d, sigma_i, threshold, 
        strategy, cells, Nselect, precision, verbose
      );

  if(result)
  {
    int nx, ny, nz; 
    
    printf("Read image\n");
  
    float *Ic=iio_read_image_float_vec(image, &nx, &ny, &nz);

    if(verbose)
      printf(
        "Parameters:\n"
        "  input image: '%s', output image: '%s', output corner file: %s\n"
        "  gaussian: %d, gradient: %d, measure: %d, K: %f, sigma_d: %f  \n"
        "  sigma_i: %f, threshold: %f, strategy: %d, Number of cells: %d\n"
        "  Nselect: %d, precision: %d, nx: %d, ny: %d, nz: %d\n",
        image, out_image, out_file, gaussian, gradient, measure, k, 
        sigma_d, sigma_i, threshold, strategy, cells, Nselect, 
        precision, nx, ny, nz
      );

    if (Ic!=NULL)
    {
      std::vector<harris_corner> corners;
      float *I=new float[nx*ny];
      
      printf("Convert image to greyscale levels\n");

      if(nz>1)
        rgb2gray(Ic, I, nx, ny, nz);
      else
        for(int i=0;i<nx*ny;i++)
          I[i]=Ic[i];

      struct timeval start, end;

      if (verbose)
        gettimeofday(&start, NULL);


      printf("Execute Harris method\n");

      //compute Harris' corners
      harris(
        I, corners, gaussian, gradient, measure, k, sigma_d, 
        sigma_i, threshold, strategy, cells, Nselect, 
        precision, nx, ny, verbose
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
        printf("Draw points\n");

        draw_points(Ic, corners, strategy, cells, nx, ny, nz, 2*sigma_i);
        
        printf("Save image\n");

        iio_save_image_float_vec(out_image, Ic, nx, ny, nz);
      }

      if(out_file!=NULL)
      {
        FILE *fd=fopen(out_file,"w");
        
        printf("Write points to file\n");

       
        fprintf(fd, "Number of points: %ld\n", corners.size());
        for(unsigned int i=0;i<corners.size();i++)
          fprintf(fd, "%f %f %f\n", corners[i].x, corners[i].y, corners[i].Mc);
        fclose(fd);
      }

      printf("Free memory\n");

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
