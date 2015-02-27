//   SpinOp.h -*- C++ -*- Operators for a spin or two-level system.
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

#ifndef _SpinOp_hhh
#define _SpinOp_hhh 1

#include "PrimOp.h"

class SigmaX: public PrimaryOperator
{
public:
  SigmaX() : PrimaryOperator(0,SPIN) {};
  SigmaX(int freedom) : PrimaryOperator(freedom,SPIN) {};
  virtual void applyTo(State&,int,double);
};

class SigmaY: public PrimaryOperator
{
public:
  SigmaY() : PrimaryOperator(0,SPIN) {};
  SigmaY(int freedom) : PrimaryOperator(freedom,SPIN) {};
  virtual void applyTo(State&,int,double);
};

class SigmaZ: public PrimaryOperator
{
public:
  SigmaZ() : PrimaryOperator(0,SPIN) {};
  SigmaZ(int freedom) : PrimaryOperator(freedom,SPIN) {};
  virtual void applyTo(State&,int,double);
};

class SigmaPlus: public PrimaryOperator
{
public:
  SigmaPlus() : PrimaryOperator(0,SPIN) {};
  SigmaPlus(int freedom) : PrimaryOperator(freedom,SPIN) {};
  virtual void applyTo(State&,int,double);
};

class SigmaMinus: public PrimaryOperator
{
public:
  SigmaMinus(int = 0) {
    error("SigmaMinus not implemented. Use SigmaPlus::hc().");
  }
  virtual void applyTo(State&,int,double) {};
};

#endif
