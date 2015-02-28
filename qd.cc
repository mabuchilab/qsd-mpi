// Test program for qd modules.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <fftw3.h>

#include "Coeff.h"
#include "Constants.h"
#include "Mesh.h"

extern "C" {
    void fcoeff_(int* n, double* c);
}
main()
{
  Coeff cc;
  Mesh  mm;
  double* c;
  int i;
  double* rho;
  double* lapl;
  double* exl;
  double* dif;
  double alpha = 3.0, r2;

//  fftw_complex in[2], out[2];
	fftw_complex* in = (fftw_complex*)fftw_malloc(2*sizeof(fftw_complex));
	fftw_complex* out = (fftw_complex*)fftw_malloc(2*sizeof(fftw_complex));

  fftw_plan p;
  p = fftw_plan_dft_1d(2, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);

exit(1);

  i = cc.getN();
//  c = new  double[2*cc.getN()+1];
    c = (double*)malloc( (2*i+1)*sizeof(double) );
   fcoeff_(&i, c);
   r2 = 0.0; 
   for (i=0;i<=2*cc.getN();i++) {
    printf("%d, %12.6f\n", i, c[i]);
    }

//  Test Laplacian module
  int n = mm.getN();
  rho = (double*)malloc( n*n*sizeof(double) );
  lapl = (double*)malloc( n*n*sizeof(double) );
  exl = (double*)malloc( n*n*sizeof(double) );
  dif = (double*)malloc( n*n*sizeof(double) );


  int ix, iy;
  for (ix=0;ix<n;ix++) {
    for (iy=0;iy<n;iy++) {
       r2 = (-1)*(pow(mm.x(ix,iy), 2.0) + pow(mm.y(ix,iy), 2.0)/(alpha*alpha));
       rho[ix*n+iy] = exp(r2)/(pi*alpha*2.0);
       printf("%d  %d  r2  %12.6f rho %20.6f %d ndiv2  %6.4f\n", n, ix*n+iy, r2, rho[ix*n+iy], n, (n/2.0)*mm.getDelta() );
    }
  }
//  printf("rho is loaded\n");
  mm.laplacian(rho, lapl, c, cc.getN());


  // Exact values of Laplacian
   for (ix=0;ix<n;ix++) {
    for (iy=0;iy<n;iy++) {
       r2 = (pow(mm.x(ix,iy), 2.0) + pow(mm.y(ix,iy), 2.0));
       exl[ix*n+iy] = (4.0/pow(alpha,2.0))*((r2/pow(alpha,2.0))-1.0)*rho[(ix*n)+iy];
//       exl[ix*n+iy] = (-1*2.0/pow(alpha,2.0))*(mm.x(ix,iy)+mm.y(ix,iy))*rho[(ix*n)+iy]; // first order derivative
       dif[ix*n+iy] = lapl[ix*n+iy] - exl[ix*n+iy];
       printf("inx %d  diff  %20.6f \n", ix*n+iy, dif[ix*n+iy]);
    }
  }
   printf("Total error %20.6f \n", mm.dotproduct(dif, dif));

//  ~Coeff();
}


