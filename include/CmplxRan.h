//   CmplxRan.h -*- C++ -*- Complex random number generators.
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

#ifndef _CmplxRan_hhh
#define _CmplxRan_hhh

#include <cstdlib>
#include "Complex.h"
#include "RNG.h"

class ComplexRandom
{
public:
  virtual ~ComplexRandom() {};
  virtual Complex operator()() = 0;
};

class ComplexZero : public ComplexRandom
{
public:
  ComplexZero() {};
  virtual Complex operator()() { return Complex(0,0); };
};

class ComplexUniform : public ComplexRandom
{
public:
  ComplexUniform(RNG *gen);
  virtual Complex operator()();
private:
  RNG *pGenerator;
};

class ComplexNormal : public ComplexRandom
{
public:
  ComplexNormal(RNG *gen);
  virtual Complex operator()();
private:
  RNG *pGenerator;
};

class ComplexNormalTest : public ComplexRandom
{ // Simple portable random number generator for test purposes only.
public:
  ComplexNormalTest(long seed);
  virtual Complex operator()();
private:
  long idum;
};

class ComplexQuadratic : public ComplexRandom
{
public:
  ComplexQuadratic(RNG *gen);
  virtual Complex operator()();
private:
  RNG *pGenerator;
  double sigma;        // sigma^2 = variance
  double xMax;         // P(x) = pMax - x^2  (-xMax < x < xMax)
  double pMax;         //  (not normalized)
};

class RefinedComplexNormal : public ComplexRandom
{
public:
  RefinedComplexNormal(RNG *gen1, RNG *gen2, int theOffset=1);
  virtual ~RefinedComplexNormal();
  virtual Complex operator()();
private:
  RNG *generator1;
  RNG *generator2;
  int offset;
  int cacheSize;
  int cacheCount;
  Complex* cache;
};

#endif
