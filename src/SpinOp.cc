//   SpinOp.cc -- Operators for a spin or two-level system.
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

#include "SpinOp.h"

static const char rcsid[] = "$Id: SpinOp.cc,v 3.1 1996/11/19 10:05:08 rschack Exp $";

void SigmaX::applyTo(State& v, int hc, double)
//
// `hc=HC' and `hc=NO_HC' are identical for Hermitian operators.
{
  Complex tmp = v[1];
  v[1] = v[0];
  v[0] = tmp;
}

void SigmaY::applyTo(State& v, int hc, double)
//
// `hc=HC' and `hc=NO_HC' are identical for Hermitian operators.
{
  Complex tmp = v[1];
  v[1] = v[0];
  v[1].timesMinusI();
  v[0] = tmp;
  v[0].timesI();
}

void SigmaZ::applyTo(State& v, int hc, double)
//
// `hc=HC' and `hc=NO_HC' are identical for Hermitian operators.
{
  v[0] *= -1;
}

void SigmaPlus::applyTo(State& v, int hc, double)
{
  switch( hc ) {
  case NO_HC:
    v[1] = v[0];
    v[0] = 0;
    break;
  case HC:
    v[0] = v[1];
    v[1] = 0;
    break;
  default:
    error("Unknown option in SigmaPlus::applyTo.");
  }
}
