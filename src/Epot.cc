//   Mesh.cc -- Grid defined on Real space.
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <math.h>
#include <stdlib.h>

#include "Complex.h"
#include "Mesh.h"
#include "Epot.h"

// These are global variables of file scope, used to avoid recalculating
// the square roots of the integers more often than necessary.
double omega;

Epot::Epot() {    //  Constructor for uninitialized Mesh

 omega = 0.22;
}

Epot::~Epot() {
}

double* Epot::harmonic()
{
  double* v;
  int ix, iy;
  double r2;

  Mesh mm;
  int n = mm.getN();
  v = new double[2*n];
  for (ix=0;ix<n;ix++) {
    for (iy=0;iy<n;iy++) {
       r2 =  pow(mm.x(ix, iy), 2.0) + pow(mm.y(ix, iy), 2.0);
       v[ix*n+iy] = 0.5*pow(omega, 2.0)*r2;
    }
  }
  return v;
}

double* Epot::quartic()
{
  double* v;
  int ix, iy;
  double r2, alpha;

  Mesh mm;
  int n = mm.getN();
  v = new double[2*n];
  alpha = 0.00008;
  for (ix=0;ix<n;ix++) {
    for (iy=0;iy<n;iy++) {
       r2 =  pow(mm.x(ix, iy), 2.0) + pow(mm.y(ix, iy), 2.0);
       v[ix*n+iy] = alpha*r2*r2;  //power 4
    }
  }
  return v;
}

