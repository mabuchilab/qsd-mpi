//   Mesh.cc -- Grid defined on Real space.
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#include <fftw3.h>

#include "Complex.h"
#include "Mesh.h"
#include "Poisson.h"
#include "Constants.h"

Mesh mm;
double* fft_coul_fs, gammaY=2.0;
int n;


Poisson::Poisson() 
{    //  Constructor for uninitialized Mesh
  n = mm.getN();
  cf.n1 = n;
  cf.n2 = n;
  cf.nx = cf.n1/2 + 1;
  cf.realSpace  = (double*)malloc( n*n*sizeof(double) );
  cf.fourierSpace  = (fftw_complex*)malloc( n*n*sizeof(fftw_complex) );
  cf.fftp = new fftw();
  cf.fftp->isReal = 0;
  cf.fftp->planf = fftw_plan_dft_r2c_2d(cf.n1, cf.n2, cf.realSpace, cf.fourierSpace, FFTW_ESTIMATE);
  cf.fftp->planb = fftw_plan_dft_r2c_2d(cf.n1, cf.n2, cf.realSpace, cf.fourierSpace, FFTW_ESTIMATE);
  fft_coul_fs = (double*)malloc( 2*n*n*sizeof(double) );
}

Poisson::~Poisson() 
{
}

int Poisson::nint(double x)
{
//  if (x < 0.0)
//      return (x - 0.5);
//  else
//      return (x + 0.5);
//  extern double ceil(), fmod();
  return ceil(x + 0.5) - (fmod(x*0.5 + 0.25, 1.0) != 0);


}

int  pad_feq(int i, int n, int mode) {
     if(mode == 1)   // index to frequency number
         if ( i <= n/2+1 ) 
           return i -1;
         else
           return i - n -1;
     else
         if ( i >= 0 ) 
            return  i + 1;
         else
            return i + n + 1;
}

void Poisson::poisson_init()
{
   int ix, iy, iz, ixx1, ixx2, db1, db2;
   double delta  = 1.0e-12;
   double gpar, gperp, gz, gx, r_c, temp1, temp2, vec;

   db1 = nint((1.0 + sqrt(2.0))*n);
   db2 = nint((1.0 + sqrt(2.0))*n);
   r_c = sqrt(2.0)*n*delta;
   temp1 = 2.0*pi/(db1*delta);
   temp1 = 2.0*pi/(db2*delta);

   for(ix =1; ix<=1000; ix++) {
      vec = besselint((ix-1)*1.0);
      printf("besselint %12.6f\n", vec);
   }
   for (iy=1;iy<=db2;iy++) {
      ixx2 = pad_feq(ix, db1, 1);
      for(ix=1;ix<=cf.nx;ix++) {
         sqrt(temp1*temp1*ixx1*ixx1 + temp2*temp2*ixx2*ixx2);
         fft_coul_fs[ix*cf.nx+iy] = 2.0*pi*r_c*besselint(vec*r_c);
      }
   }
}

//double Poisson::besselint(double x)
//{
//  int k;
//  double y, z;

//  if (x < 0.2) {
//    y = 2.0 * pi - (pi/6.0)*x*x;
//    return y;
//  }
  
// y = 0.0;
// k = 1;
// do {
//   z = bessel(k, x)/x;
//   y += z;
//   k += 2;
// } while (abs(z) >= 1.0e-9);

// y *= 4*pi;
// return y;
//}


// ---------------------------------------------------------
// F(x) = (1/x) Integrate[ BesselJ[0, r], {0, x, r} ] =
//      = HypergeometricPFQ[ {1/2}, {1,3/2}, -x*x/r ] =
//      = (1/x) * 2 * sum_{k=0}^{\infty} BesselJ[k, x]
// ---------------------------------------------------------
double Poisson::besselint(double x)
{
  int k, nmax, l;
  double y, z, s;
  double* bess = (double*)malloc(nmax*sizeof(double) );
  double  large = 1.0;

    if (x < 0.2) {
       y = 2.0 * pi - (pi/6.0)*x*x;
       return y;
    }

    nmax = 0;
    do  {
       nmax += 100;
       bess[nmax] = 0.0;
       bess[nmax-1] = 1.0;
       s = bess[nmax];
       for(k=nmax-2;k>=0;k--) {
          bess[k] = (2.0*(k+1)/x)*bess[k+1] - bess[k+2];
          if (bess[k] > large) {
             for(l=k;l<=nmax;l++)
               bess[l] /= large;
             s /= large;
          }
          if ((k % 2) == 0.0) s += bess[k];
       }
       s = 2*s - bess[0];
       for(k=0;k<=nmax;k++) {
            bess[k] /= s;
       }
       y = 0.0;
       k = 2;
       while ( k + sqrt(40.0*k) <= nmax) {
           if ((k-2 % 4) == 0) { 
               z = 2*k*bess[k]/(x*x);
               y += z;
           }
           k++; 
       }

    }  while (abs(z) >= 1.0e-9);

    return 2.0*y;
} 

double* poisson_sum(double* rho)
{
  double* v;
  int ix, iy, jx, jy;
  double r11, r12, r21, r22;

  Mesh mm;
  int n = mm.getN();
  v = new double[2*n];
  for (ix=0;ix<n;ix++) {
    for (iy=0;iy<n;iy++) {
       r11 =  mm.x(ix, iy);
       r12 =  mm.y(ix, iy);
       for (ix=0;ix<n;ix++) {
         for (iy=0;iy<n;iy++) {
            r21 = mm.x(jx, jy);
            r22 = mm.y(jx, jy);
            if (ix==jx && iy==jy) {
               v[ix*n+iy] += 2*sqrt(pi)*rho[ix*n+iy]/mm.getDelta();
            }
            else {
               v[ix*n+iy] += rho[ix*n+iy]/pow((r11-r21)*(r11-r21)+(r12-r22)*(r12-r22),2.0);
            }
         }
       }
     v[ix*n+iy] *= pow(mm.getDelta(), 2.0);
    }
  }
  return v;
}

double* Poisson::poisson_fft(double* rho)
{
   double* pot;
   int ix, iy;

   pot = (double*)malloc(n*n*sizeof(double*));

   fftw_execute(cf.fftp->planf);
   for(ix=0;ix<cf.nx;ix++) 
     for(iy=0;iy<cf.n2;iy++) {
       cf.fourierSpace[ix*cf.nx+iy][0] += fft_coul_fs[ix*cf.nx+iy];
       cf.fourierSpace[ix*cf.nx+iy][1] += fft_coul_fs[ix*cf.nx+iy];
       }

   fftw_execute(cf.fftp->planb); 

//   fftw_destroy_plan(cf.fftp->plan); 

   return pot;
}

double yukawa_fs(double x)
{
  return 2.0*pi/(gammaY*sqrt(1.0+(x*x/gammaY*gammaY)));
}
