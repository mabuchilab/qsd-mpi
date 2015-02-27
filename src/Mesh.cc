//   Mesh.cc -- Grid defined on Real space.
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <math.h>
#include <stdlib.h>

#include "Complex.h"
#include "Mesh.h"


// These are global variables of file scope, used to avoid recalculating
// the square roots of the integers more often than necessary.

 const int n=81;  // total number of points in each direction

 double mlen; // mesh length
 double delta;
 double temp;

Mesh::Mesh() {    //  Constructor for uninitialized Mesh

  mlen =  50.0;  // the length L of the cube orig len 50.8
  delta = mlen/(double)(n-1);   // spacing between the points

}

Mesh::~Mesh() {
}

double Mesh::x(int ix, int iy)
{
  double  delx = ((ix -1) - (n/2.0))*delta;
  return delx;
}

double Mesh::y(int ix, int iy)
{
  double dely = ((iy -1) - (n/2.0))*delta;
  return dely;
}

double Mesh::dotproduct(double* x1, double* x2)
{
  int i;
  double sum;

  sum = 0.0;
  for (i=0;i<n*n;i++) {
    sum += x1[i]*x2[i];
  }
  printf("%20.6f   dot product\n", sum);
  sum *= delta*delta;
  return sum;
}

int Mesh::getN() 
{

  return n;

}

double Mesh::getDelta()
{
   
   return delta;

}

void Mesh::laplacian(double* f, double* lapl, double *c, int i)
{
  int ix, iy, k;


  for (ix=0;ix<n;ix++) {
    for (iy=0;iy<n;iy++) {
      lapl[(ix*n)+iy] = 0.0;
        lapl[ix*n+iy]  += c[0]*f[ix*n+iy];
        for (k=1;k<=i;k++) {
           if(iy+k>=1 && iy+k<n) lapl[ix*n+iy]  += c[2*i+1-k]*f[ix*n+iy+k];
           if(ix+k>=1 && ix+k<n) lapl[ix*n+iy]  += c[2*i+1-k]*f[(ix+k)*n+iy];
           if(iy+k>=1 && iy+k<n) lapl[ix*n+iy]  += c[i+1-k]*f[ix*n+iy+k];
           if(ix+k>=1 && ix+k<n) lapl[ix*n+iy]  += c[i+1-k]*f[(ix+k)*n+iy];
        }
        lapl[ix*n+iy] /= delta*delta;
    }
  }
}

void zlaplacian(Complex** f, Complex** lapl, double *c)
{
  int ix, iy, k;


  for (ix=1;ix<=n;ix++) {
    for (iy=1;iy<=n;iy++) {
      lapl[ix][iy] = temp;
        for (k=-4;k<=4;k++) {
           if(iy+k>=1 && iy+k<=n) lapl[ix][iy]  += c[k]*f[ix][iy+k];
           if(ix+k>=1 && ix+k<=n) lapl[ix][iy]  += c[k]*f[ix+k][iy];
        }
      }
    }
  for (ix=1;ix<=n;ix++) {
    for (iy=1;iy<=n;iy++) {
       lapl[ix][iy] /= pow(delta,2.0);
    }
  }
}
