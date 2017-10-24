#ifndef _HARRIS_H
#define _HARRIS_H

#include<vector>

void harris(
  float *I,
  std::vector<float> &x,
  std::vector<float> &y,
  float alpha,
  float sigma_i,
  float sigma_n,
  int   radius,
  float percentage,
  int   nobel_measure,
  int   precision,
  int   nx,
  int   ny,
  int   verbose,
  int   forensics
);

#endif


/*
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
*/
