//   Coeff.cc -- This object provides methods to compute the coefficients
//   for the second derivative in a regular on-dimensional grid (Laplacian)
/////////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "Coeff.h"


extern "C" {
    void fcoeff_(int* n, double* c);
}

Coeff::Coeff()  {    //  Constructor for uninitialized Mesh
  m = 2;  // the order of the derivative
  del = 1.0;
  n = 40;  // approximation for second order derivative
  int i;
  x = new double[2*n+1];
  // it assigns x(1)=delta,x(2)=2delta...x(5)=-delta...x(8)=-4*delta
  for (i=0;i<=n;i++) {
    x[i] = i*del;  
  }
  for (i=n+1;i<=2*n;i++) {
    x[i] = (n-i)*del;
  }
}

//Coeff::~Coeff() { free(); }             // Destructor

int Coeff::getN() {

  return n;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// FUNCTION COEFF
// ==============

// INPUT:
//   m [integer] : the order of the derivatives representation.
//   n [integer] : The number of points given is 2*n
//   x [real(8), dimension(2*n)] : positions of the points. The "problem" point position is not given,
//     and assumed to be zero.
// ---------
// OUTPUT:
//   c [real(8), dimension(2*n+1)] : the coefficients of the points. The first one corresponds to the
//     the coefficient at the problem points (which is always minus the sum of all the others), whereas
//     the rest are ordered in the same manner that were given in array x.
//   coeff [integer] : error code. It is the error code of the LAPACK subroutine dgels

// Calculates the coefficients for the representation of the m-th order
// derivative of a function at a given point, given that we will have access
// to the values of this function at 2*n points around it (besides the value
// of this function at the problem point). Typically this means n points to the
// left, and n points to the right, but that is not mandatory.

// NOTES:
// ------
// It requires BLAS dgesv subroutine.
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void Coeff::coefficients(int n, double *c)
{
  fcoeff_(&n, c);
}

void Coeff::coefficients(double *c)
{
  int i, j, k, lwork, info, mm, nn, lda, ldb, nrhs;
  char* trans;
  double* a;
  double* e;
  double* work;
  double wkopt;
  double sum;

  a = new double[2*n*2*n]; 
  e = new double[2*n];
  
  trans = new char[12];
  trans = "No transpose";
  nn = 2*n;
  mm = 2*n;
  lda = 2*n;
  ldb = 2*n;
  nrhs = 1;

  for (i=0;i<2*n;i++) {
    e[i] = 0.0;
    for (j=0;j<2*n;j++) {
      a[i+j] = pow(x[j+1],(double)(i+1));
//      printf("%d   %12.6f\n", i+j, a[i+j]);
      }
    }
  k = 1;
  for (i=1;i<=2*n;i++) {
      k = k*i;
      if (m==i) {
       e[i-1] = k;
//       break;
      }
       printf("%d %d   %12.6f\n", m, i, e[i-1]);
  }

  lwork = -1;
  printf("%d  %d  %d  %d %d %d %d Before calling dgels_1\n", mm,nn,nrhs,lda,ldb,lwork,info);
//  dgels(trans, &mm, &nn, &nrhs, a, &lda, e, &ldb, &wkopt, &lwork, &info);
//  printf("after calling dgela_1 lwork\n"); //  %12.6f  %d  %12.6f\n", lwork, info, wkopt);
//  lwork = (int)wkopt;
//  work = (double*)malloc( lwork*sizeof(double) );
//  dgels_("No transpose", &mm, &nn, &nrhs, a, &nn, e, &oo, work, &lwork, &info);

  
//  for (i=1;i<2*n;i++) {
//   c[i] = e[j];
//   sum += e[j];
//  }
//  c[0] = -1*sum;
  printf("Assigned all cs\n");
  /* Free workspace */
//  free( (void*)work );
}

