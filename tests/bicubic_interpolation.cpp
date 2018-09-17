// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2018, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
// All rights reserved.


#include "bicubic_interpolation.h"

#include<cmath>

//types of transformations
#define TRANSLATION_TRANSFORM 2
#define EUCLIDEAN_TRANSFORM   3
#define SIMILARITY_TRANSFORM  4
#define AFFINITY_TRANSFORM    6
#define HOMOGRAPHY_TRANSFORM  8



/**
 *
 *  Function to invert a transforms 
 *  W(x;p) = W(x;p1)^-1
 *
 */
void inverse_transform
(
  float *p1,        //input transform
  float *p,         //output transform
  const int nparams  //number of parameters
)
{
  switch(nparams) 
  {
    default: case TRANSLATION_TRANSFORM: //p=(tx, ty)
      p[0]=-p1[0];
      p[1]=-p1[1];
      break;
    case EUCLIDEAN_TRANSFORM:  //p=(tx, ty, tita)
    {
      const float a=p1[0];
      const float b=p1[1];
      const float c=p1[2];
      p[0]=-a*cos(c)-b*sin(c);
      p[1]= a*sin(c)-b*cos(c);
      p[2]=-c;
      break;
    }
    case SIMILARITY_TRANSFORM:  //p=(tx, ty, a, b)
    {
      const float a=p1[2];
      const float b=p1[3];
      const float c=p1[0];
      const float d=p1[1];
      const float det=2*a+a*a+b*b+1;
      if(det*det>1E-10)
      {
        p[0]=(-c-a*c-b*d)/det;
        p[1]=(-d-a*d+b*c)/det;
        p[2]=(a+1)/det-1;
        p[3]=-b/det;
      }
      else p[0]=p[1]=p[2]=p[3]=0;
    }
    break;
    case AFFINITY_TRANSFORM:    //p=(tx, ty, a00, a01, a10, a11)
    {
      const float a=p1[2];      
      const float b=p1[3];
      const float c=p1[0];
      const float d=p1[4];
      const float e=p1[5];
      const float f=p1[1];
      const float det=a-b*d+e+a*e+1;    
      if(det*det>1E-10)
      {
        p[0]=(-c+b*f-c*e)/det;
        p[1]=(-f-a*f+c*d)/det;
        p[2]=(e+1)/det-1;
        p[3]=-b/det;
        p[4]=-d/det;
        p[5]=(a+1)/det-1; 
      }
      else p[0]=p[1]=p[2]=p[3]=p[4]=p[5]=0;
    }
    break;
    case HOMOGRAPHY_TRANSFORM:   //p=(h00, h01,..., h21)
    {
      const float a=p1[0];
      const float b=p1[1];
      const float c=p1[2];
      const float d=p1[3];
      const float e=p1[4];
      const float f=p1[5];
      const float g=p1[6];
      const float h=p1[7];
      const float det=(-a+b*d-e-a*e-1);
      if(det*det>1E-10)
      {
        p[0]=(f*h-e-1)/det-1;
        p[1]=(b-c*h)/det;
        p[2]=(c-b*f+c*e)/det;
        p[3]=(d-f*g)/det;
        p[4]=(-a+c*g-1)/det-1;
        p[5]=(f+a*f-c*d)/det;
        p[6]=(g-d*h+g*e)/det;
        p[7]=(h+a*h-b*g)/det;
      }
      else p[0]=p[1]=p[2]=p[3]=p[4]=p[5]=p[6]=p[7]=0;
    }
    break;

  }
}


/**
 *
 *  Function to transform a 2D point (x,y) through a parametric model
 *
 */
void project
(
  float x,   //x component of the 2D point
  float y,   //y component of the 2D point
  float *p,  //parameters of the transformation
  float &xp, //x component of the transformed point
  float &yp, //y component of the transformed point
  int nparams //number of parameters
)
{
  switch(nparams) {
    default: case TRANSLATION_TRANSFORM: //p=(tx, ty) 
      xp=x+p[0];
      yp=y+p[1];
      break;
    case EUCLIDEAN_TRANSFORM:   //p=(tx, ty, tita)
      xp=cos(p[2])*x-sin(p[2])*y+p[0];
      yp=sin(p[2])*x+cos(p[2])*y+p[1];
      break;
    case SIMILARITY_TRANSFORM:  //p=(tx, ty, a, b)
      xp=(1+p[2])*x-p[3]*y+p[0];
      yp=p[3]*x+(1+p[2])*y+p[1];
      break;
    case AFFINITY_TRANSFORM:    //p=(tx, ty, a00, a01, a10, a11)
      xp=(1+p[2])*x+p[3]*y+p[0];
      yp=p[4]*x+(1+p[5])*y+p[1];
      break;
    case HOMOGRAPHY_TRANSFORM:  //p=(h00, h01,..., h21)
      double d=p[6]*x+p[7]*y+1;
      xp=((1+p[0])*x+p[1]*y+p[2])/d;
      yp=(p[3]*x+(1+p[4])*y+p[5])/d;
      break;
  }
}



/**
  *
  * Neumann boundary condition test
  *
**/
int
neumann_bc (int x, int nx, bool & out)
{
  if (x<0)
    {
      out = true;
      x = 0;
    }
  else if (x >= nx)
    {
      out = true;
      x = nx - 1;
    }
  return x;
}


/**
  *
  * Bicubic interpolation in one dimension
  *
**/
float
cubic_interpolation(
  float v[4],  //interpolation points
  float x      //point to be interpolated
)
{
  return v[1] + 0.5 * x * (v[2] - v[0]
                           + x * (2.0 * v[0] - 5.0 * v[1] + 4.0 * v[2] - v[3]
                                  + x * (3.0 * (v[1] - v[2]) + v[3] - v[0])));
}


/**
  *
  * Bicubic interpolation in two dimension
  *
**/
float
bicubic_interpolation(
  float p[4][4], //array containing the interpolation points
  float x,       //x position to be interpolated
  float y        //y position to be interpolated
)
{
  float v[4];
  v[0] = cubic_interpolation (p[0], y);
  v[1] = cubic_interpolation (p[1], y);
  v[2] = cubic_interpolation (p[2], y);
  v[3] = cubic_interpolation (p[3], y);
  return cubic_interpolation (v, x);
}


/**
  *
  * Compute the bicubic interpolation of a point in an image. 
  * Detects if the point goes outside the image domain
  *
**/
float
bicubic_interpolation(
  float *input,//image to be interpolated
  float uu,    //x component of the vector field
  float vv,    //y component of the vector field
  int nx,       //width of the image
  int ny,       //height of the image
  bool border_out //if true, put zeros outside the region
)
{
  int sx = (uu < 0) ? -1 : 1;
  int sy = (vv < 0) ? -1 : 1;

  int x, y, mx, my, dx, dy, ddx, ddy;
  bool out = false;

  x = neumann_bc ((int) uu, nx, out);
  y = neumann_bc ((int) vv, ny, out);
  mx = neumann_bc ((int) uu - sx, nx, out);
  my = neumann_bc ((int) vv - sy, ny, out);
  dx = neumann_bc ((int) uu + sx, nx, out);
  dy = neumann_bc ((int) vv + sy, ny, out);
  ddx = neumann_bc ((int) uu + 2 * sx, nx, out);
  ddy = neumann_bc ((int) vv + 2 * sy, ny, out);

  if (out && border_out) 
    return 0;
  else
    {
      //obtain the interpolation points of the image
      float p11 = input[mx  + nx * my];
      float p12 = input[x   + nx * my];
      float p13 = input[dx  + nx * my];
      float p14 = input[ddx + nx * my];

      float p21 = input[mx  + nx * y];
      float p22 = input[x   + nx * y];
      float p23 = input[dx  + nx * y];
      float p24 = input[ddx + nx * y];

      float p31 = input[mx  + nx * dy];
      float p32 = input[x   + nx * dy];
      float p33 = input[dx  + nx * dy];
      float p34 = input[ddx + nx * dy];

      float p41 = input[mx  + nx * ddy];
      float p42 = input[x   + nx * ddy];
      float p43 = input[dx  + nx * ddy];
      float p44 = input[ddx + nx * ddy];

      //create array
      float pol[4][4] = { 
        {p11, p21, p31, p41}, {p12, p22, p32, p42},
        {p13, p23, p33, p43}, {p14, p24, p34, p44}
      };

      //return interpolation
      return bicubic_interpolation (pol, (float) uu - x, (float) vv - y);
    }
}





/**
  *
  * Compute the bicubic interpolation of an image from a parametric trasform
  *
**/
void bicubic_interpolation(
  float *input,   //image to be warped
  std::vector<int> &p, //selected points
  float *output,  //warped output image with bicubic interpolation
  float *params,  //x component of the vector field
  int nparams,     //number of parameters of the transform
  int nx,          //width of the image
  int ny,          //height of the image 
  bool border_out  //if true, put zeros outside the region
)
{
  for (unsigned int i=0; i<p.size(); i++)
  {
    float x, y;
    float x1=p[i]%nx;
    float y1=(int)(p[i]/nx);

    //transform coordinates using the parametric model
    project(x1, y1, params, x, y, nparams);
    
    //obtain the bicubic interpolation at position (uu, vv)
    output[i]=bicubic_interpolation(input, x, y, nx, ny, border_out);
  }
}


/**
  *
  * Compute the bicubic interpolation of an image from a parametric trasform
  *
**/
void bicubic_interpolation(
  float *input,   //image to be warped
  float *output,  //warped output image with bicubic interpolation
  float *params,  //x component of the vector field
  int nparams,     //number of parameters of the transform
  int nx,          //width of the image
  int ny,          //height of the image 
  bool border_out  //if true, put zeros outside the region
)
{
  for (int i=0; i<ny; i++)
    for (int j=0; j<nx; j++)
    {
      float x, y;

      //transform coordinates using the parametric model
      project(j, i, params, x, y, nparams);
      
      //obtain the bicubic interpolation at position (uu, vv)
      output[i*nx+j]=bicubic_interpolation(
        input, x, y, nx, ny, border_out
      );
    }
}

