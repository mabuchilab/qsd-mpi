//   CmplxRan.cc -- Complex Random number generators.
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

#include <math.h>
#include "CmplxRan.h"

#ifndef M_SQRT1_2
#define M_SQRT1_2   0.70710678118654752440
#endif

ComplexUniform::ComplexUniform(RNG *gen)
{
  pGenerator = gen;
}

Complex ComplexUniform::operator()()
{
  return Complex(pGenerator->asDouble(), 0.0 );
}

ComplexNormal::ComplexNormal(RNG *gen)
{
  pGenerator = gen;
}

Complex ComplexNormal::operator()()
{
  for(;;) {
    double u1 = pGenerator -> asDouble();
    double u2 = pGenerator -> asDouble();
    double v1 = 2 * u1 - 1;
    double v2 = 2 * u2 - 1;
    double w = (v1 * v1) + (v2 * v2);
	    
    if (w <= 1) {
      double y = sqrt( (-2 * log(w)) / w);
      return M_SQRT1_2 * Complex(v1*y, v2*y);  
                                // M_SQRT1_2 = 1/sqrt(2) [cf math.h]
    }
  }
}

ComplexNormalTest::ComplexNormalTest(long seed)
{
  idum = seed;
}

Complex ComplexNormalTest::operator()()
{  // Numerical Recipes in C, Second Edition, p. 279.
  const long IA=16807;
  const long IM=2147483647;
  const double AM=(1.0/IM);
  const long IQ=127773;
  const long IR=2836;
  const long MASK=123459876;

  idum ^= MASK;

  for(;;) {
    long k = idum/IQ;
    idum = IA*(idum-k*IQ)-IR*k;
    if (idum < 0) idum += IM;
    double u1 = AM * idum;

    k = idum/IQ;
    idum = IA*(idum-k*IQ)-IR*k;
    if (idum < 0) idum += IM;
    double u2 = AM * idum;

    double v1 = 2 * u1 - 1;
    double v2 = 2 * u2 - 1;
    double w = (v1 * v1) + (v2 * v2);
	    
    if (w <= 1) {
      double y = sqrt( (-2 * log(w)) / w);
      return M_SQRT1_2 * Complex(v1*y, v2*y);  
                                // M_SQRT1_2 = 1/sqrt(2) [cf math.h]
    }
  }
  idum ^= MASK;
}

RefinedComplexNormal::RefinedComplexNormal(RNG *g1,RNG *g2,int theOffset)
{
  if( theOffset < 1 ) {
    cerr << "Error in RefinedComplexNormal: non-positive offset." << endl;
    exit(1);
  }
  offset = theOffset;
  cacheSize = 2*offset;
  cache = new Complex[cacheSize];
  cacheCount = cacheSize;
  generator1 = g1;
  generator2 = g2;
}

RefinedComplexNormal::~RefinedComplexNormal()
{
#ifndef NON_GNU_DELETE
  delete[] cache;
#else
  delete[cacheSize] cache;
#endif
}

Complex RefinedComplexNormal::operator()()
{
  if (cacheCount < cacheSize)
    return cache[cacheCount++];
  else {
    Complex* xi1;
    Complex* xi2;
    xi1 = new Complex[offset];
    xi2 = new Complex[offset];
    int i;
    for( i=0; i<offset; i++ ) {
      for(;;) {
	double u1 = generator1 -> asDouble();
	double u2 = generator1 -> asDouble();
	double v1 = 2 * u1 - 1;
	double v2 = 2 * u2 - 1;
	double w = (v1 * v1) + (v2 * v2);
	if (w <= 1) {
	  double y = sqrt( (-2 * log(w)) / w);
	  Complex xi(v1*y, v2*y);
	  xi1[i] = xi;
	  break;
	}
      }
      for(;;) {
	double u1 = generator2 -> asDouble();
	double u2 = generator2 -> asDouble();
	double v1 = 2 * u1 - 1;
	double v2 = 2 * u2 - 1;
	double w = (v1 * v1) + (v2 * v2);
	if (w <= 1) {
	  double y = sqrt( (-2 * log(w)) / w);
	  Complex xi(v1*y, v2*y);
	  xi2[i] = xi;
	  break;
	}
      }
    }
    for( i=0; i<offset; i++ ) {
      cache[i] = 0.5 * (xi1[i]+xi2[i]);
      cache[offset+i] = 0.5 * (xi1[i]-xi2[i]);
    }
    cacheCount = 0;
    return cache[cacheCount++];
  }
}

ComplexQuadratic::ComplexQuadratic(RNG *gen)
{
  pGenerator = gen;
  sigma = 1.0/sqrt(2);
  xMax = sqrt(5) * sigma;
  pMax = xMax*xMax;        // P(x) = pMax - x*x  (-xMax < x < xMax)
}

Complex ComplexQuadratic::operator()()
{
  double c[2];
  for( int i=0; i<2; i++ ) 
  for(;;) {
    double u1 = pGenerator -> asDouble();
    double u2 = pGenerator -> asDouble();
    double x = xMax*( 2*u1 - 1 );
    double p = pMax * u2;
    if ( p  <=  pMax - x*x ) {
      c[i] = x;
      break;
    }
  }
  return Complex( c[0],c[1] );
}
