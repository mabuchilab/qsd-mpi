//   FieldOp.cc -- Operators for a harmonic oscillator mode.
//
//   Copyright (C) 1995  Todd Brun and Ruediger Schack
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
//   ----------------------------------------------------------------------
//   If you improve the code or make additions to it, or if you have
//   comments or suggestions, please contact us:
//
//   Dr. Todd Brun			        Tel    +44 (0)171 775 3292
//   Department of Physics                      FAX    +44 (0)181 981 9465
//   Queen Mary and Westfield College           email  t.brun@qmw.ac.uk
//   Mile End Road, London E1 4NS, UK
//
//   Dr. Ruediger Schack                        Tel    +44 (0)1784 443097
//   Department of Mathematics                  FAX    +44 (0)1784 430766
//   Royal Holloway, University of London       email  r.schack@rhbnc.ac.uk
//   Egham, Surrey TW20 0EX, UK
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib>

#include "FieldOp.h"

static const char rcsid[] = "$Id: FieldOp.cc,v 3.1 1996/11/19 10:05:08 rschack Exp $";

static int intSqrtSize = 0;
static double* intSqrt = 0;
static double* halfIntSqrt = 0;

static void computeIntSqrt(int newSize)
  //
  // Compute integer square roots.
{
  if(intSqrtSize < newSize) {
#ifndef NON_GNU_DELETE
    delete[] intSqrt;
    delete[] halfIntSqrt;
#else
    delete[intSqrtSize] intSqrt;
    delete[intSqrtSize] halfIntSqrt;
#endif
    intSqrtSize = newSize;
    intSqrt = new double[newSize];
    halfIntSqrt = new double[newSize];
    for( int i=0; i<newSize; i++ ) {
      intSqrt[i] = sqrt(i);
      halfIntSqrt[i] = sqrt(0.5*i);
    }
  }
}

void FieldTransitionOperator::applyTo(State& v, int hc, double)
//  |i><j|  (Transition from |j> to |i>)
{
#ifdef DEBUG_TRACE
  std::cout << "FieldTransitionOperator(freedom="<<getFreedom()<<", i="<<i<<", j="<<j<<")::applyTo(v, hc="<<hc<<") entered." << std::endl;
#endif
  int vSize = v.size();  // v.size() is the number of elements of the vector v.
#ifndef OPTIMIZE_QSD
  if ( i >= vSize || j >= vSize || i < 0 || j < 0 )
    error("Transition out of range in FieldTransitionOperator::applyTo.");
#endif
  int k;
  Complex center = v.centerCoords();
  if( center != 0 ) {
    error("FieldTransitionOperator is incompatible with moving basis.");
  }
  else {
    switch( hc ) {
    case NO_HC:
      v[i] = v[j];
      for( k=0; k<i; k++ )
        v[k] = 0;
      for( k=i+1; k<vSize; k++ )
        v[k] = 0;
      break;
    case HC:
      v[j] = v[i];
      for( k=0; k<j; k++ )
        v[k] = 0;
      for( k=j+1; k<vSize; k++ )
        v[k] = 0;
      break;
    default:
      error("Unknown option in FieldTransitionOperator::applyTo.");
    }
  }
}

void AnnihilationOperator::applyTo(State& v, int hc, double)
{
#ifdef DEBUG_TRACE
  std::cout << "AnnihilationOperator(freedom="<<getFreedom()<<")::applyTo(v, hc="<<hc<<") entered." << std::endl;
#endif
  int i;
  Complex center = v.centerCoords();
  int vSize = v.size();
  if( intSqrtSize < vSize )
    computeIntSqrt(vSize);
  if( center != 0 ) {
    //    error("Moving basis not implemented for AnnihilationOperator.");
    switch( hc ) {
      case NO_HC:
        for( i=1; i<vSize; i++ )
          v[i-1] = center*v[i-1] + intSqrt[i] * v[i];
        v[vSize-1] = center*v[vSize-1];
        break;
      case HC:
        center = conj(center);
        for( i=vSize-1; i>0; i-- )
          v[i] = center*v[i] + intSqrt[i] * v[i-1];
        v[0] = center*v[0];
        break;
      default:
        error("Unknown option AnnihilationOperator::applyTo.");
    }
  }
  else {
    switch( hc ) {
      case NO_HC:
        for( i=1; i<vSize; i++ )
          v[i-1] = intSqrt[i] * v[i];
        v[vSize-1] = 0;
        break;
      case HC:
        for( i=vSize-1; i>0; i-- )
          v[i] = intSqrt[i] * v[i-1];
        v[0] = 0;
        break;
      default:
        error("Unknown option AnnihilationOperator::applyTo.");
    }
  }
}

void LocalLower::applyTo(State& v, int hc, double)
{
#ifdef DEBUG_TRACE
  std::cout << "LocalLower(freedom="<<getFreedom()<<")::applyTo(v, hc="<<hc<<") entered." << std::endl;
#endif
  int i;
  int vSize = v.size();
  if( intSqrtSize < vSize )
    computeIntSqrt(vSize);
  switch( hc ) {
    case NO_HC:
      for( i=1; i<vSize; i++ )
        v[i-1] = intSqrt[i] * v[i];
      v[vSize-1] = 0;
      break;
    case HC:
      for( i=vSize-1; i>0; i-- )
        v[i] = intSqrt[i] * v[i-1];
      v[0] = 0;
      break;
    default:
      error("Unknown option LocalLower::applyTo.");
  }
}

void NumberOperator::applyTo(State& v, int hc, double)
  //
  // `hc=HC' and `hc=NO_HC' are identical for Hermitian operators.
{
#ifdef DEBUG_TRACE
  std::cout << "NumberOperator(freedom="<<getFreedom()<<")::applyTo(v, hc="<<hc<<") entered." << std::endl;
#endif
  int i;
  Complex center = v.centerCoords();
  int vSize = v.size();
  if( intSqrtSize < vSize )
    computeIntSqrt(vSize);
  if( center != 0 ) {
    //    error("Moving basis not implemented for AnnihilationOperator.");
    Complex previous,temp,vzero,vnth;
    Complex cstar=conj(center);
    double mag=norm(center);
    vzero = mag*v[0] + cstar*v[1];
    vnth = (vSize+mag-1)*v[vSize-1] + center*intSqrt[vSize-1]*v[vSize-2];
    previous = v[0];
    for (i=1; i<vSize-1; i++) {
      temp = v[i];
      v[i] = (i+mag)*v[i] + cstar*intSqrt[i+1]*v[i+1] + center*intSqrt[i]*previous;
      previous = temp;
    }
    v[0] = vzero;
    v[vSize-1] = vnth;
  }
  else {
    for( i=0; i<vSize; i++ )
      v[i] *= i;
  }
}

void XOperator::applyTo(State& v, int hc, double)
  //
  // `hc=HC' and `hc=NO_HC' are identical for Hermitian operators.
{
#ifdef DEBUG_TRACE
  std::cout << "XOperator(freedom="<<getFreedom()<<")::applyTo(v, hc="<<hc<<") entered." << std::endl;
#endif
  Complex center = v.centerCoords();
  int vSize = v.size();
  Complex previous,temp,vzero,vnth;
  if( intSqrtSize < vSize )
    computeIntSqrt(vSize);
  if( center != 0 ) {
    double xcent = intSqrt[2]*real(center);
    vzero = xcent*v[0] + halfIntSqrt[1]*v[1];
    vnth = xcent*v[vSize-1] + halfIntSqrt[vSize-1]*v[vSize-2];
    previous = v[0];
    for (int i=1; i<vSize-1; i++) {
      temp = v[i];
      v[i] = xcent*v[i] + halfIntSqrt[i+1]*v[i+1] + halfIntSqrt[i]*previous;
      previous = temp;
    }
    v[0] = vzero;
    v[vSize-1] = vnth;
  }
  else {
    vzero = halfIntSqrt[1]*v[1];
    vnth =  halfIntSqrt[vSize-1]*v[vSize-2];
    previous = v[0];
    for (int i=1; i<vSize-1; i++) {
      temp = v[i];
      v[i] = halfIntSqrt[i+1]*v[i+1] + halfIntSqrt[i]*previous;
      previous = temp;
    }
    v[0] = vzero;
    v[vSize-1] = vnth;
  }
}

void POperator::applyTo(State& v, int hc, double)
  //
  // `hc=HC' and `hc=NO_HC' are identical for Hermitian operators.
{
#ifdef DEBUG_TRACE
  std::cout << "POperator(freedom="<<getFreedom()<<")::applyTo(v, hc="<<hc<<") entered." << std::endl;
#endif
  Complex center = v.centerCoords();
  int vSize = v.size();
  Complex previous,temp,vzero,vnth;
  Complex imaginary_unit(0,1);
  if( intSqrtSize < vSize )
    computeIntSqrt(vSize);
  if( center != 0 ) {
    double pcent = intSqrt[2]*imag(center);
    vzero = pcent*v[0] - imaginary_unit*halfIntSqrt[1]*v[1];
    vnth = pcent*v[vSize-1] + imaginary_unit*halfIntSqrt[vSize-1]*v[vSize-2];
    previous = v[0];
    for (int i=1; i<vSize-1; i++) {
      temp = v[i];
      v[i] = halfIntSqrt[i]*previous - halfIntSqrt[i+1]*v[i+1];
      v[i].timesI();
      v[i] += pcent*temp;
      previous = temp;
    }
    v[0] = vzero;
    v[vSize-1] = vnth;
  }
  else {
    vzero = - imaginary_unit*halfIntSqrt[1]*v[1];
    vnth = imaginary_unit*halfIntSqrt[vSize-1]*v[vSize-2];
    previous = v[0];
    for (int i=1; i<vSize-1; i++) {
      temp = v[i];
      v[i] = halfIntSqrt[i]*previous - halfIntSqrt[i+1]*v[i+1];
      v[i].timesI();
      previous = temp;
    }
    v[0] = vzero;
    v[vSize-1] = vnth;
  }
}

void DisplacementOperator::applyTo(State& v, int hc, double)
{
#ifdef DEBUG_TRACE
  std::cout << "DisplacementOperator(freedom="<<getFreedom()<<", alpha="<<alpha<<")::applyTo(v, hc="<<hc<<") entered." << std::endl;
#endif
  int i,j,k,n,m;
  Complex center = v.centerCoords();
  int vSize = v.size();
  if( intSqrtSize < vSize )
    computeIntSqrt(vSize);

  if( matrixSize < vSize ) {

    for( k=0; k<matrixSize; k++ )       // reallocation of matrix
#ifndef NON_GNU_DELETE
      delete[] matrix[k];
#else
    delete[matrixSize] matrix[k];
#endif
#ifndef NON_GNU_DELETE
    delete[] matrix;
    delete[] vv;
#else
    delete[matrixSize] matrix;
    delete[matrixSize] vv;
#endif
    matrixSize = vSize;
    vv = new Complex[vSize];
    matrix = new Complex*[vSize];
    for( k=0; k<matrixSize; k++ )
      matrix[k] = new Complex[vSize];

    double nalpha = norm(alpha);        // Computation of matrix
    double c = exp( -0.5*nalpha );
    for( n=0; n<matrixSize; n++ ) {
      for( m=0; m<=n; m++ ) {
        Complex s=c;
        for( k=0; k<n-m; k++ )
          s *= alpha*intSqrt[n-k]/(k+1);
        Complex sum=s;
        for( k=1; k<=m; k++ ) {
          s *= -nalpha*(m-k+1) / (k*(k+n-m));
          sum += s;
        }
        matrix[n][m] = sum;
      }
    }
    for( n=0; n<matrixSize; n++ ) {
      for( m=n+1; m<matrixSize; m++ ) {
        if( 2*((m-n)/2) == m-n )
          matrix[n][m] = conj(matrix[m][n]);
        else
          matrix[n][m] = -conj(matrix[m][n]);
      }
    }
  }
  switch( hc ) {
    case NO_HC:
      for( i=0; i<vSize; i++ ) {
        vv[i] = 0;
        for( j=0; j<vSize; j++ )
          vv[i] += matrix[i][j] * v[j];
      }
      for( i=0; i<vSize; i++ )
        v[i] = vv[i];
      break;
    case HC:
      for( i=0; i<vSize; i++ ) {
        vv[i] = 0;
        for( j=0; j<vSize; j++ )
          vv[i] += conj(matrix[j][i]) * v[j];
      }
      for( i=0; i<vSize; i++ )
        v[i] = vv[i];
      break;
    default:
      error("Unknown option DisplacementOperator::applyTo.");
  }
  if( center != 0 ) {
    double phi=2*imag(conj(center)*alpha);
    Complex phase( cos(phi), sin(phi) );
    if( hc == HC ) phase = conj(phase);
    for( i=0; i<vSize; i++ )
      v[i] *= phase;
  }
}

void RealDisplacementOperator::applyTo(State& v, int hc, double)
{
#ifdef DEBUG_TRACE
  std::cout << "RealDisplacementOperator(freedom="<<getFreedom()<<", alpha="<<alpha<<")::applyTo(v, hc="<<hc<<") entered." << std::endl;
#endif
  int i,j,k,n,m;
  Complex center = v.centerCoords();
  int vSize = v.size();
  if( intSqrtSize < vSize )
    computeIntSqrt(vSize);

  if( matrixSize < vSize ) {

    for( k=0; k<matrixSize; k++ )       // reallocation of matrix
#ifndef NON_GNU_DELETE
      delete[] matrix[k];
#else
    delete[matrixSize] matrix[k];
#endif
#ifndef NON_GNU_DELETE
    delete[] matrix;
    delete[] vv;
#else
    delete[matrixSize] matrix;
    delete[matrixSize] vv;
#endif
    matrixSize = vSize;
    vv = new Complex[vSize];
    matrix = new double*[vSize];
    for( k=0; k<matrixSize; k++ )
      matrix[k] = new double[vSize];

    double nalpha = alpha*alpha;        // Computation of matrix
    double c = exp( -0.5*nalpha );
    for( n=0; n<matrixSize; n++ ) {
      for( m=0; m<=n; m++ ) {
        double s=c;
        for( k=0; k<n-m; k++ )
          s *= alpha*intSqrt[n-k]/(k+1);
        double sum=s;
        for( k=1; k<=m; k++ ) {
          s *= -nalpha*(m-k+1) / (k*(k+n-m));
          sum += s;
        }
        matrix[n][m] = sum;
      }
    }
    for( n=0; n<matrixSize; n++ ) {
      for( m=n+1; m<matrixSize; m++ ) {
        if( 2*((m-n)/2) == m-n )
          matrix[n][m] = matrix[m][n];
        else
          matrix[n][m] = -matrix[m][n];
      }
    }
  }
  switch( hc ) {
    case NO_HC:
      for( i=0; i<vSize; i++ ) {
        vv[i] = 0;
        for( j=0; j<vSize; j++ )
          vv[i] += matrix[i][j] * v[j];
      }
      for( i=0; i<vSize; i++ )
        v[i] = vv[i];
      break;
    case HC:
      for( i=0; i<vSize; i++ ) {
        vv[i] = 0;
        for( j=0; j<vSize; j++ )
          vv[i] += matrix[j][i] * v[j];
      }
      for( i=0; i<vSize; i++ )
        v[i] = vv[i];
      break;
    default:
      error("Unknown option RealDisplacementOperator::applyTo.");
  }
  if( center != 0 ) {
    double phi=2*imag(conj(center)*alpha);
    Complex phase( cos(phi), sin(phi) );
    if( hc == HC ) phase = conj(phase);
    for( i=0; i<vSize; i++ )
      v[i] *= phase;
  }
}
