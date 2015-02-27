//   AtomOp.cc -- Operators for an atom or N-level system
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

#include "AtomOp.h"

static const char rcsid[] = "$Id: AtomOp.cc,v 3.1 1996/11/19 10:05:08 rschack Exp $";

void TransitionOperator::applyTo(State& v, int hc, double)
//  |i><j|  (Transition from |j> to |i>)
{
  int vSize = v.size();  // v.size() is the number of elements of the vector v.
#ifndef OPTIMIZE_QSD
  if ( i >= vSize || j >= vSize || i < 0 || j < 0 )
    error("Transition out of range in TransitionOperator::applyTo.");
#endif
  int k;
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
    error("Unknown option in TransitionOperator::applyTo.");
  }
}
