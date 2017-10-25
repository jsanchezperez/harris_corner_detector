#include "interpolation.h"

#define MAX_ITERATIONS 100


#include <stdio.h>

/**
  *
  * Evaluate the function in (x,y)
  * 
  * f(x)=a0 x^2y^2 + a1 x^2y + a2 xy^2 + a3 x^2 + a4 y^2 +
  *      a5 xy + a6 x + a7 y + a8
  *
**/
float f(
  float *a,
  float x,
  float y
)
{
  return a[0]*x*x*y*y+a[1]*x*x*y+a[2]*x*y*y+
         a[3]*x*x+a[4]*y*y+a[5]*x*y+a[6]*x+a[7]*y+a[8];
}


/**
  *
  * Compute the coefficients of the interpolation function
  * 
  * f(x)=a0 x^2y^2 + a1 x^2y + a2 xy^2 + a3 x^2 + a4 y^2 +
  *      a5 xy + a6 x + a7 y + a8
  *
**/
void polynomial_coefficients(
  float *M,
  float *a
)
{
  a[0]=M[4]-0.5*(M[1]+M[3]+M[5]+M[7])+0.25*(M[0]+M[2]+M[6]+M[8]);
  a[1]=0.5*(M[1]-M[7])+0.25*(-M[0]-M[2]+M[6]+M[8]);
  a[2]=0.5*(M[3]-M[5])+0.25*(-M[0]+M[2]-M[6]+M[8]);
  a[3]=0.5*(M[3]+M[5])-M[4];
  a[4]=0.5*(M[1]+M[7])-M[4];
  a[5]=0.25*(M[0]-M[2]-M[6]+M[8]);
  a[6]=0.5*(M[5]-M[3]);
  a[7]=0.5*(M[7]-M[1]);
  a[8]=M[4];
}


/**
  *
  * Compute the Gradient of the interpolation function
  * 
  *
**/
void polynomial_gradient(
  float dx, //current x-position
  float dy, //current y-position
  float *a, //polynomial coefficients
  float *D  //output Hessian
)
{
  D[0]=2*a[0]*dx*dy*dy+2*a[1]*dx*dy+2*a[2]*dy*dy+2*a[3]*dx+a[5]*dy+a[6];
  D[1]=2*a[0]*dx*dx*dy+2*a[1]*dx*dx+2*a[2]*dx*dy+2*a[4]*dy+a[5]*dx+a[7];
}


/**
  *
  * Compute the Hessian of the interpolation function
  * 
  *
**/
void Hessian(
  float dx, //current x-position
  float dy, //current y-position
  float *a, //polynomial coefficients
  float *H  //output Hessian
)
{
  H[0]=2*a[0]*dy*dy+2*a[1]*dy+2*a[3];
  H[1]=4*a[0]*dx*dy+2*a[1]*dx+2*a[2]*dy+a[5];
  H[2]=2*a[0]*dx*dx+2*a[2]*dx+2*a[4];
}


/**
  *
  * Solve the system H^-1*D
  * 
  *
**/
bool solve(
  float *H, //Hessian
  float *D, //gradient
  float *b  //output increment
)
{
  float det=H[0]*H[2]-H[1]*H[1];
  
  if(det*det<10E-10)
    return false;
  else{
    b[0]=(D[0]*H[2]-D[1]*H[1])/det;
    b[1]=(D[1]*H[0]-D[0]*H[1])/det;
    return true;
  }
}


/**
  *
  * Apply Newton method to find maximum of the interpolation function
  *
**/
bool maximum_interpolation(
  float *M, //values of the surfare (9 values)
  float &x, //x solution
  float &y, //y solution
  float TOL //stopping criterion threshold
)
{
  float D[2],b[2],H[3],a[9];
  
  float dx=0, dy=0;
  polynomial_coefficients(M,a);
  
  int i=0;
  do
  {
    //compute the gradient and Hessian in the current position
    polynomial_gradient(dx,dy,a,D);
    Hessian(dx,dy,a,H);
    
    //solve the system for estimating the next increment
    if(!solve(H,D,b)) 
      return false;
    
    //move the current position
    dx-=b[0];
    dy-=b[1];        
    
    i++;
  }
  while(D[0]*D[0]+D[1]*D[1]>TOL && i<MAX_ITERATIONS);
    
  //printf("Iterations> %d  increment: (%f, %f)\n", i, dx, dy);
  
  //check that (dx, dy) are inside the allowed boundaries 
  if(dx>1 || dx<-1 || dy>1 || dy<-1)
    return false;
  else
  {
    x+=dx; y+=dy;
    return true;
  }
}

float G(float x, float y)
{
 float mx=-0.95;
 float my=0.95;
 float r=(-(x-mx)*(x-mx)-(y-my)*(y-my));
 printf("f(%f, %f)=%f\n", x, y, r);
 return r;
}


/*void main()
{
 float M[9],x=0,y=0;
 int i;
 M[0]=G(-1.,-1.) ;
 M[1]=G(0.,-1.);
 M[2]=G(1.,-1.);

 M[3]=G(-1.,0.);
 M[4]=G(0.,0.);
 M[5]=G(1.,0.);

 M[6]=G(-1.,1.);
 M[7]=G(0.,1.);
 M[8]=G(1.,1.);

 if(maximum_interpolation(M, x, y))
     printf("Resultado = %f, %f \n",x,y);

}*/

