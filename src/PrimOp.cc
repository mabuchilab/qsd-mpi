//   PrimOp.cc -- Base class for special operators.
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

#include <stdlib.h>
#include <iostream.h>

#include "PrimOp.h"

static const char rcsid[] = "$Id: PrimOp.cc,v 3.1 1996/11/19 10:05:08 rschack Exp $";

int PrimaryOperator::pflag = 0;       // Flag set in constructor.

void PrimaryOperator::resetFreedom(int freedom)
//
// Reset the freedom on which the PrimaryOperator acts
{
#ifndef OPTIMIZE_QSD
  if( freedom < 0 ) {
    cerr << "Illegal freedom passed to PrimaryOperator::resetFreedom!" << endl;
    exit(1);
  }
#endif
  myFreedom = freedom;
}

void IdentityOperator::applyTo(State&,int,double)
{
#ifdef DEBUG_TRACE
  cout << "IdentityOperator::applyTo entered." << endl;
#endif
}

void NullOperator::applyTo(State& v,int,double)
{
#ifdef DEBUG_TRACE
  cout << "NullOperator::applyTo entered." << endl;
#endif
  for( int i=0; i<v.size(); i++ )
    v[i] = 0;
}

/////////////////////////////////////////////////////////////////////////////
// The following `ExampleOperator' can be used as a template to create new //
// special operators derived from the `PrimaryOperator' class.             //
/////////////////////////////////////////////////////////////////////////////

void ExampleOperator::applyTo(State& v, int hc, double t)
{
  // `applyTo' defines how `ExampleOperator' operates on a vector `v'. On 
  // exit, `v' must contain the result. The vector components are accessed
  // using square brackets as in `v[i]'.
  //
  // `ExampleOperator' does not use a moving basis. For examples using a
  // moving basis, see the include file `FieldOp.h'.

  Complex im(0,1);
  int vSize = v.size();  // v.size() is the number of elements of the vector v.

  double tp = t*parameter;   // For time-independent operators, ignore t.
  int k;

  switch( hc ) {    // This switch is needed only for non-Hermitian operators.
                    // For Hermitian operators, ignore 'hc'.

  case NO_HC:                      // Here, define how ExampleOperator operates
    for( k=0; k<vSize; k++ )       // on the vector `v'.
      v[k] = exp(im*tp*k) * v[k];
    break;

  case HC:                         // Here, define how the Hermitian conjugate
    for( k=0; k<vSize; k++ )       // of ExampleOperator operates on `v'.
      v[k] = exp(-im*tp*k) * v[k];
    break;

  default:
    error("Unknown option in ExampleOperator::applyTo.");
  }
}

/////////// End of template /////////////////////////////////////////////////

