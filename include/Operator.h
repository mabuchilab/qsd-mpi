//   Operator.h -*- C++ -*- Operator algebra in Hilbert space.
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

#ifndef _Operator_hhh
#define _Operator_hhh 1

#include "State.h"
#include "Complex.h"

// The `Operator' class and its derived classes use the following definitions
// and functions of the `State' class:
//
// enum FreedomType;
// enum ImaginaryUnit;
// const ImaginaryUnit IM;
// const ImaginaryUnit M_IM;
// double size();                       // Return number of vector components.
// void xerox(State psi);               // Allocate memory and copy vector
//                                      // components to an unitialized state.
// State& operator*=(const Complex&);   // Multiplication by a complex scalar.
// State& operator*=(double);           // Multiplication by a real scalar.
// State& operator*=(ImaginaryUnit);    // Multiplication by IM or M_IM.
// State& operator+=(const State&);     // Addition
// State& operator-=(const State&);     // Subtraction

typedef double (*RealFunction) (double); 
typedef Complex (*ComplexFunction) (double); 
//
// This defines the types `RealFunction' and `ComplexFunction' 
// which can be used as in
// `RealFunction sinus = sin; double y = sinus(2.0);'
// or in
// `Complex si(double t) { return sin(t); } 
//  ComplexFunction sinus = si; Complex y = sinus(2.0);'

enum Command{ OPERATOR, OPERATOR_HC, COMPLEX, REAL, IMAG, M_IMAG,
	      CFUNC, CCFUNC, RFUNC, PLUS, MINUS, TIMES, UNINITIALIZED };

class Operator
{
public:

  Operator();                           // Default constructor.
  Operator(const Operator&);            // Copy constructor.
  virtual ~Operator();                  // Destructor.
  Operator& operator=(const Operator&); // Assignment.

  State operator*(const State& psi) const;
  friend State& operator*=(State& psi, const Operator& X);
  //
  // Application to a state.
  // The `*=' form is more efficient. In `psi1 = X * psi;', a temporary
  // `State' object is created in addition to `psi' and `psi1'. In
  // `psi *= X;', however, no additional `State' object is created.

  Operator operator+(const Operator&) const;
  Operator operator-(const Operator&) const;
  Operator operator*(const Operator&) const;
  Operator& operator+=(const Operator&);
  Operator& operator-=(const Operator&);
  // Addition, subtraction, and multiplication of two operators.
  // `*=' not implemented because of ambiguity.

  friend Operator operator-(const Operator&);
  // The negative of an operator. `X = -Y;' is equivalent to `X = (-1)*Y;'.
  friend Operator operator+(const Operator&);

  Operator pow(int n) const;   // Integer power of operator (n > 0).

  Operator hc() const;         // Hermitian conjugate.

  Operator operator*(const Complex&) const;
  Operator operator*(double) const;
  Operator operator*(ImaginaryUnit) const;
  Operator& operator*=(const Complex&);
  Operator& operator*=(double);
  Operator& operator*=(ImaginaryUnit);
  friend Operator operator*(const Complex&,const Operator&);
  friend Operator operator*(double,const Operator&);
  friend Operator operator*(ImaginaryUnit,const Operator&);
  // Scalar multiplication. Multiplication by a real number and multiplication
  // by the imaginary unit `im' are encoded separately for efficiency. The
  // constant `im' and the type `ImaginaryUnit' are defined in "State.h".
  // The constant `im' can only be used to directly multiply a `State' or an
  // `Operator', i.e., expressions like `Complex c = 3*im;' are illegal.

  void printCommandStack() const;  // For debugging only.

  //-------- Time-dependent operators ------------

  Operator& operator()(double t); 
  // Specify `time=t' when applying a time-dependent `Operator' to a `State'.
  // `time' is passed to each time-dependent `PrimaryOperator', 
  // `ComplexFunction', and `RealFunction' that is part of the time-dependent
  // `Operator'.
  // Calling a time-dependent `Operator' without the `time' argument should
  // be avoided.
  // If a time-dependent `Operator' is applied without giving the `time' 
  // argument, the same time is used as in the last call. 
  // The default is `time=0'. Example:
  //
  // double t = 5.0;
  // IdentityOperator Id;
  // State psi(10);   
  // RealFunction f = cos;
  // Operator Cos = Id*f;
  // State psi1 = Cos * psi;    // Multiply psi by cos(0)=1.
  // State psi1 = Cos(t)*psi;   // Multiply psi by cos(t).
  // State psi1 = Cos * psi;    // Multiply psi by cos(t).

  Operator operator*(ComplexFunction) const;
  Operator operator*(RealFunction) const;
  Operator& operator*=(ComplexFunction);
  Operator& operator*=(RealFunction);
  friend Operator operator*(ComplexFunction,const Operator&);
  friend Operator operator*(RealFunction,const Operator&);
  // Multiplication by a `ComplexFunction' or a `RealFunction'.

protected:
  
  enum HcSwitch { HC, NO_HC };
  // Constants used to indicate if a primary
  // operator or its Hermitian conjugate is applied.
  void error( const char* ) const;

           #ifndef NON_GNU_INHERIT
private:
           #endif

  friend class PrimaryOperator;    // This gives `PrimaryOperator' access to 
                                   // the private data and functions.
  // Private type definitions:
  struct Stack {
    Command* com;          // The command stack.
    PrimaryOperator** op;  // Stack of pointers to PrimaryOperator
    Complex* c;            // Stack of complex scalars.
    double* r;             // Stack of real scalars.
    ComplexFunction* cf;   // Stack of complex-valued functions.
    RealFunction* rf;      // Stack of real-valued functions.
    int comSize;           // Size of command stack.
    int opSize;            // Size of PrimaryOperator stack.
    int cSize;             // Size of complex scalar stack.
    int rSize;             // Size of real scalar stack.
    int cfSize;            // Size of complex function stack.
    int rfSize;            // Size of real function stack.
  };
  struct StackPtr {
    int com;
    int op;
    int c;
    int r;
    int cf;
    int rf;
  };

  // Private data:
  Stack stack;
  double time;
  static int flag;    // Set by `Operator' constructors, tested by 
                      // `PrimaryOperator' constructors. 

  // Private functions:
  void eval( StackPtr&, State& ) const;
  void offsetCopyStack( Operator&, int, int, int, int, int, int ) const;
  void copy( const Operator& );
  void allocate( int,int,int,int,int,int );
  void deallocate();
  Operator dagger( StackPtr& ) const;
  Operator(char*);
};

#endif
