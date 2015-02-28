//   PrimOp.h -*- C++ -*- Base class for special operators.
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

// Contains the classes 
// PrimaryOperator, IdentityOperator, NullOperator, ExampleOperator

#ifndef _PrimOp_hhh
#define _PrimOp_hhh 1

#include "Operator.h"

// PrimaryOperator is an abstract class, which means that no instances
// of PrimaryOperator can be created. PrimaryOperator serves as an
// interface to "special operators" derived from the PrimaryOperator class. 
// The classes IdentityOperator, NullOperator, and ExampleOperator defined
// below are examples of special operators.
//
class PrimaryOperator: public Operator 
{
public:

  PrimaryOperator() {};   // Default constructor.

  PrimaryOperator(int freedom, FreedomType type) : Operator("") {
  //
  // This constructor is called by the constructors of any special operator
  // derived from PrimaryOperator. 
  // PrimaryOperator can set the stack because PrimaryOperator
  // is declared a friend by its base class `Operator'.
  //
    if( flag && !pflag ) {
      pflag = 1;
      std::cerr << "Warning: It is recommended to define all special operators,\n"
	<< "i.e., all instances of classes derived from `PrimaryOperator',\n"
        << "e.g., instances of `IdentityOperator' or `NullOperator',\n"
	<< "before an instance of the `Operator' class is defined.\n"
        << "This avoids the following potential problem:\n"
        << "Any instance of an `Operator' contains pointers to one or more\n"
        << "special operators. If one of these special operators goes out of\n"
        << "scope before the instance of the `Operator' class itself, then\n"
        << "the latter points to an unallocated memory location, which can\n"
        << "be catastrophic. A typical consequence is a segmentation fault.\n";
    }
    stack.com[0] = OPERATOR;
    stack.op[0] = this;         // Points to itself.
    myFreedom = freedom;
    myType = type;
  };

  virtual void applyTo(State&, int, double) = 0;
  // Abstract virtual function. Prototypes `applyTo' for any special operator
  // derived from `PrimaryOperator'. 
  // In addition, precludes creation of instances of `PrimaryOperator'.

  void resetFreedom(int);
  // Resets the freedom number on which the primary operator acts

private:                  // The private data are inherited by any special 
                          // operator derived from `PrimaryOperator'.

  int myFreedom;        
  FreedomType myType;
  friend class Operator;  // Gives `Operator' access to myFreedom and myType.
  static int pflag;
};

class IdentityOperator: public PrimaryOperator 
{
public:
  IdentityOperator() : PrimaryOperator(0,ALL) {};
  IdentityOperator(int freedom) : PrimaryOperator(freedom,ALL) {};
  virtual void applyTo(State&,int,double);
};

class NullOperator: public PrimaryOperator 
{
public:
  NullOperator() : PrimaryOperator(0,ALL) {};
  NullOperator(int freedom) : PrimaryOperator(freedom,ALL) {};
  virtual void applyTo(State&,int,double);
};

/////////////////////////////////////////////////////////////////////////////
// The following `ExampleOperator' can be used as a template to create new //
// special operators derived from the `PrimaryOperator' class.             //
/////////////////////////////////////////////////////////////////////////////

class ExampleOperator: public PrimaryOperator
{
public:

  ExampleOperator() : PrimaryOperator(0,ALL) { parameter = 0; };
  //
  // Default constructor. ExampleOperator operates on a degree of freedom
  // of type ALL. Other possible types are enumerated under `FreedomType'
  // in the include file `State.h'.
  // The default constructor should always use `0' as the first argument in
  // the call to the `PrimaryOperator' constructor, since operators are
  // expected to operate on the first degree of freedom by default.
  // ExampleOperator has a private parameter `parameter' which is
  // initialized to `0' by default.

  ExampleOperator(double p,int freedom=0) : PrimaryOperator(freedom,ALL) {
    parameter = p;
  };
  // Constructor that explicitely sets `parameter' and `myFreedom'.

  virtual void applyTo(State& v, int hc, double t);
  // Member function which contains the actual code defining ExampleOperator.

private:         // There could be more private data, or none at all.

  double parameter;
};

////////// End of template. /////////////////////////////////////////////////

#endif
