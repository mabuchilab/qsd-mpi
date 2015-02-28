//   Operator.cc -- Operator algebra in Hilbert space
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

#include <cstdlib>
#include <iostream>

#include "Operator.h"
#include "PrimOp.h"

static const char rcsid[] = "$Id: Operator.cc,v 3.1 1996/11/19 10:05:08 rschack Exp $";

int Operator::flag = 0;   // Static flag, set by constructors, tested by
                          // PrimaryOperator class.

void Operator::error( const char* msg ) const
{
  std::cerr << "Fatal error in class Operator or derived class:\n  " << msg << std::endl;
  exit(1);
}

void Operator::allocate
  (int comSize,int opSize=0,int cSize=0,int rSize=0,int cfSize=0,int rfSize=0)
//
// Private function, used in object construction and assignment.
{
  if( comSize<0 || opSize<0 || cSize<0 || rSize<0 || cfSize<0 || rfSize<0 )
    error("Negative stack size in allocate.");
  stack.com = 0;
  stack.op = 0;
  stack.c = 0;
  stack.r = 0;
  stack.cf = 0;
  stack.rf = 0;
  if( (stack.comSize = comSize) > 0) stack.com = new Command[comSize];
  if( (stack.opSize = opSize) > 0) stack.op = new PrimaryOperator*[opSize];
  if( (stack.cSize = cSize) > 0) stack.c = new Complex[cSize];
  if( (stack.rSize = rSize) > 0) stack.r = new double[rSize];
  if( (stack.cfSize = cfSize) > 0) stack.cf = new ComplexFunction[cfSize];
  if( (stack.rfSize = rfSize) > 0) stack.rf = new RealFunction[rfSize];
}

void Operator::deallocate()
//
// Private function, used in object destruction and assignment.
{
  #ifndef NON_GNU_DELETE
    delete[] stack.com;
    delete[] stack.op;
    delete[] stack.c;
    delete[] stack.r;
    delete[] stack.cf;
    delete[] stack.rf;
  #else  //--------NON GNU CODE--------
    delete[stack.comSize] stack.com;
    delete[stack.opSize] stack.op;
    delete[stack.cSize] stack.c;
    delete[stack.rSize] stack.r;
    delete[stack.cfSize] stack.cf;
    delete[stack.rfSize] stack.rf;
  #endif //--------NON GNU CODE--------
}

void Operator::copy(const Operator& rhs)
//
// Private function, used in object construction and assignment.
{
  time = rhs.time;
  allocate( rhs.stack.comSize, rhs.stack.opSize, rhs.stack.cSize, 
	   rhs.stack.rSize, rhs.stack.cfSize, rhs.stack.rfSize );
  int i;
  for (i=0; i<stack.comSize; i++) stack.com[i] = rhs.stack.com[i];
  for (i=0; i<stack.opSize; i++) stack.op[i] = rhs.stack.op[i];
  for (i=0; i<stack.cSize; i++) stack.c[i] = rhs.stack.c[i];
  for (i=0; i<stack.rSize; i++) stack.r[i] = rhs.stack.r[i];
  for (i=0; i<stack.cfSize; i++) stack.cf[i] = rhs.stack.cf[i];
  for (i=0; i<stack.rfSize; i++) stack.rf[i] = rhs.stack.rf[i];
}

Operator::Operator()   // Default constructor.
{
  allocate(1);
  stack.com[0] = UNINITIALIZED;
  time = 0;
  flag = 1;
}

//Operator::Operator(char*)  // Private constructor called by `PrimaryOperator'.
Operator::Operator(const char*)  // Private constructor called by `PrimaryOperator'.
{                          // Does not set `flag'.               
  allocate(1,1);
  stack.com[0] = UNINITIALIZED;
  time = 0;
}

Operator::Operator(const Operator& rhs)   // Copy constructor.
{
  copy(rhs);
  flag = 1;
}

Operator::~Operator()    // Destructor.
{
  deallocate();
}

Operator& Operator::operator=(const Operator& rhs)   // Assignment.
{
  deallocate();
  copy(rhs);
  return *this;
}

void Operator::printCommandStack() const      // For debugging only.
{
  std::cout << "Operator::printCommandStack():" << std::endl;
  for( int i=0; i<stack.comSize; i++ )
    std::cout << stack.com[i] << std::endl;
}

Operator& Operator::operator()(double t)
{
  time = t;
  return *this;
}

State Operator::operator*(const State& psi) const
//
// Application of the operator (i.e., `*this') to the state `psi' using the 
// syntax `State psi0 = X * psi;' where `X' is the operator.
// Defined in terms of `*='. Creates a temporary `State' object and is
// therefore inefficient. Use `*=' wherever possible.
{
  State psi1(psi);
  return psi1 *= *this;
}

State& operator*=(State& psi, const Operator& X)
//
// Application of the operator `X' to the state `psi' using the 
// syntax `psi *= X;'. The state `psi' is modified; no extra `State'
// object is created. 
{
  #ifdef DEBUG_TRACE
    std::cout << "operator*=(State&,Operator&) entered." << std::endl;
  #endif
  #ifndef NON_GNU_SCOPE
    Operator::StackPtr stackPtr = {0,0,0,0,0,0};  // Initialize stack pointers.
  #else  //----------NON GNU CODE----
    StackPtr stackPtr = {0,0,0,0,0,0};
  #endif //----------NON GNU CODE----
  X.eval( stackPtr, psi );             // Enter the recursive `Operator::eval'.
  #ifdef DEBUG_TRACE
    std::cout << "Returning from operator*=(State&,Operator&)." << std::endl;
  #endif
  return psi;
}

void Operator::eval( StackPtr& stackPtr, State& psi ) const
/////////////////////////////////////////////////////////////////////////// 
// Private function.
// Applies the "popped" stack `stack' to the state `psi', thereby modifying
// `psi', and pops the stack by incrementing the stack pointer `stackPtr'.
// "Popped stack" means -- e.g., for the command stack `stack.com' --
// that the stack level given by the value of `stackPtr.com' is interpreted
// as the bottom level of the stack.
// 
// Example:  The assignment  `C=(A+B)*4.7;', where `A' and `B' are 
// primary operators, leads to the following contents of `C.stack':
// 
//  i    C.stack.com[i]    C.stack.op[i]     C.stack.c[i]     C.stack.r[i]
// 
//  4      OP
//  3      OP
//  2      PLUS
//  1      REAL           (pointer to A)
//  0      TIMES          (pointer to B)                          4.7
//
//  (`stack.cf' and `stack.rf' are omitted to keep the example simple.)
// 
// The statement 
// psi *= C;
// then leads to the execution of the following piece of code:
// {
//   psi *= 4.7;
//   State psi1;
//   psi1.xerox(psi);
//   psi1 *= B;
//   psi *= A;
//   psi += psi1;
// }
//////////////////////////////////////////////////////////////////////////////
{
  Command c=stack.com[stackPtr.com++];   // Take command-stack content and pop.
  switch (c) {
    
  case OPERATOR: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::eval for case OPERATOR entered." << std::endl;
    #endif
    PrimaryOperator& A = *stack.op[stackPtr.op++];
    psi.apply(A,NO_HC,A.myFreedom,A.myType,time);
    // Apply to `psi' the primary operator to which `stack.op' points.
    // Pop `stack.op'
    break;
  }

  case OPERATOR_HC: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::eval for case OPERATOR_HC entered." << std::endl;
    #endif
    PrimaryOperator& A = *stack.op[stackPtr.op++];
    psi.apply(A,HC,A.myFreedom,A.myType,time);
    // Apply to `psi' the Hermitian conjugate of the primary operator to which
    // `stack.op' points. Pop `stack.op'
    break;
  }

  case COMPLEX:                     // Multiply by the content of `stack.c'.
    psi *= stack.c[stackPtr.c++];   // Pop `stack.c'.
    break;                          // `*=' defined in class State.

  case REAL:                        // Multiply by the content of `stack.r'.
    psi *= stack.r[stackPtr.r++];   // Pop `stack.r'.
    break;                          // `*=' defined in class State.

  case IMAG:            // Multiply by the imaginary unit `IM'.
    psi *= IM;          // `*=' defined in class State.
    break;

  case M_IMAG:          // Multiply by the negative imaginary unit `M_IM'.
    psi *= M_IM;        // `*=' defined in class State.
    break;

  case CFUNC:  // Multiply by the function in `stack.cf' evaluated at time t.
    psi *= stack.cf[stackPtr.cf++](time); // Pop `stack.cf'.
    break;                                // `*=' defined in class State.

  case CCFUNC:  // Multiply by the conjugate of `stack.cf' evaluated at time t.
    psi *= conj(stack.cf[stackPtr.cf++](time)); // Pop `stack.cf'.
    break;                                      // `*=' defined in class State.

  case RFUNC:   // Multiply by the function in `stack.rf' evaluated at time t.
    psi *= stack.rf[stackPtr.rf++](time);  // Pop `stack.rf'.
    break;                                 // `*=' defined in class State.

  case PLUS: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::eval for case PLUS entered." << std::endl;
    #endif
    State psi1;                     
    psi1.xerox(psi);                // Copy psi into psi1.
    eval( stackPtr, psi1 );         // Apply (popped) stack to psi1.
    eval( stackPtr, psi );          // Apply (popped) stack to psi.
    psi += psi1;                    // Add the results.
    break;                          // `+=' defined in class State.
  }  // `psi1' goes out of scope.
    
  case MINUS: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::eval for case MINUS entered." << std::endl;
    #endif
    State psi1;                     
    psi1.xerox(psi);                // Copy psi into psi1.
    eval( stackPtr, psi1 );         // Apply (popped) stack to psi1.
    eval( stackPtr, psi );          // Apply (popped) stack to psi.
    psi -= psi1;                    // Subtract the results.
    break;                          // `-=' defined in class State.
  }  // `psi1' goes out of scope.

  case TIMES:
    #ifdef DEBUG_TRACE
    std::cout << "Operator::eval for case TIMES entered." << std::endl;
    #endif
    eval( stackPtr, psi );     // Apply (popped) stack to psi.
    eval( stackPtr, psi );     // Apply (popped) stack to (modified) psi.
    break;

  case UNINITIALIZED:    
    error("Attempt to apply an uninitialized operator to a state.");

  default:
    error("eval: Unknown object in command stack.");
  }
}

void Operator::offsetCopyStack (Operator& target,
      int n_com,int n_op=0,int n_c=0,int n_r=0,int n_cf=0,int n_rf=0) const
//
// Private function, used in `+', `-', etc. Copies `stack' into
// `target.stack', offset by integers n_com, n_op, n_c, n_r, n_cf, and n_rf.
{
  int newComSize = stack.comSize + n_com;
  int newOpSize = stack.opSize + n_op;
  int newCSize = stack.cSize + n_c;
  int newRSize = stack.rSize + n_r;
  int newCfSize = stack.cfSize + n_cf;
  int newRfSize = stack.rfSize + n_rf;
  if( target.stack.comSize < newComSize || 
     target.stack.opSize < newOpSize || 
     target.stack.cSize < newCSize ||
     target.stack.rSize < newRSize ||
     target.stack.cfSize < newCfSize ||
     target.stack.rfSize < newRfSize )
    {                        // Allocate a new stack to `target' if necessary.
      target.deallocate();
      target.allocate
	(newComSize,newOpSize,newCSize,newRSize,newCfSize,newRfSize);
    }

  int i;
  for( i=0; i<stack.comSize; i++ )          // Copy stack.com offset by n_com.
  target.stack.com[i+n_com] = stack.com[i];
  
  for( i=0; i<stack.opSize; i++ )           // Copy stack.op offset by n_op.
    target.stack.op[i+n_op] = stack.op[i];

  for( i=0; i<stack.cSize; i++ )            // Copy stack.c offset by n_c.
    target.stack.c[i+n_c] = stack.c[i];

  for( i=0; i<stack.rSize; i++ )            // Copy stack.r offset by n_r.
    target.stack.r[i+n_r] = stack.r[i];

  for( i=0; i<stack.cfSize; i++ )           // Copy stack.cf offset by n_cf.
    target.stack.cf[i+n_cf] = stack.cf[i];

  for( i=0; i<stack.rfSize; i++ )           // Copy stack.rf offset by n_rf.
    target.stack.rf[i+n_rf] = stack.rf[i];
}

Operator Operator::operator+(const Operator& X) const
////////////////////////////////////////////////////////////////////////////
// Addition of two operators.
// 
// Example: If the stacks of the operators X and Y are given by
// 
//  i    X.stack.com[i]    X.stack.op[i]     X.stack.c[i]     X.stack.r[i]
// 
//  4      OP
//  3      OP
//  2      PLUS
//  1      REAL           (pointer to A)
//  0      TIMES          (pointer to B)                          4.7
// 
// and
// 
//  i    Y.stack.com[i]    Y.stack.op[i]     Y.stack.c[i]     Y.stack.r[i]
// 
//  2      OP
//  1      REAL
//  0      TIMES          (pointer to C)                          -1.0
//  
// then the stack of the `Operator Z = X + Y;' is given by 
// 
//  i    Z.stack.com[i]    Z.stack.op[i]     Z.stack.c[i]     Z.stack.r[i]
//  
//  8      OP
//  7      OP
//  6      PLUS
//  5      REAL
//  4      TIMES
//  3      OP
//  2      REAL           (pointer to A)
//  1      TIMES          (pointer to B)                           4.7
//  0      PLUS           (pointer to C)                          -1.0
//
//  (`stack.cf' and `stack.rf' are omitted to keep the example simple.)
// 
////////////////////////////////////////////////////////////////////////////
{
  if( stack.com[0] == UNINITIALIZED || X.stack.com[0] == UNINITIALIZED )
    error("Attempt to add an uninitialized operator.");
  Operator Z;
  offsetCopyStack( Z, X.stack.comSize+1, X.stack.opSize, X.stack.cSize, 
		  X.stack.rSize, X.stack.cfSize, X.stack.rfSize );
  X.offsetCopyStack( Z, 1 );
  Z.stack.com[0] = PLUS;
  return Z;
}

Operator Operator::operator-(const Operator& X) const
//
// Subtraction of two operators. Similar to `+'.
{
  if( stack.com[0] == UNINITIALIZED || X.stack.com[0] == UNINITIALIZED )
    error("Attempt to subtract an uninitialized operator.");
  Operator Z;
  offsetCopyStack( Z, X.stack.comSize+1, X.stack.opSize, X.stack.cSize, 
		  X.stack.rSize, X.stack.cfSize, X.stack.rfSize );
  X.offsetCopyStack( Z, 1 );
  Z.stack.com[0] = MINUS;
  return Z;
}

Operator Operator::operator*(const Operator& X) const
//
// Multiplication of two operators. Similar to `+'.
{
  if( stack.com[0] == UNINITIALIZED || X.stack.com[0] == UNINITIALIZED )
    error("Attempt to multiply an uninitialized operator.");
  Operator Z;
  offsetCopyStack( Z, X.stack.comSize+1, X.stack.opSize, X.stack.cSize, 
		  X.stack.rSize, X.stack.cfSize, X.stack.rfSize );
  X.offsetCopyStack( Z, 1 );
  Z.stack.com[0] = TIMES;
  return Z;
}

Operator Operator::operator*(const Complex& alpha) const
///////////////////////////////////////////////////////////////////////////
// Scalar multiplication of an operator by a complex number.
// 
// Example: If `Complex alpha (1,1);' and the stack of the operator X 
// is given by
// 
//  i    X.stack.com[i]    X.stack.op[i]     X.stack.c[i]     X.stack.r[i]
// 
//  4      OP
//  3      OP
//  2      PLUS
//  1      COMPLEX        (pointer to A)
//  0      TIMES          (pointer to B)       3 + 2i
// 
// then the stack of the `Operator Z = X * alpha;' is given by
// 
//  i    Z.stack.com[i]    Z.stack.op[i]     Z.stack.c[i]     Z.stack.r[i]
//   
//  6      OP
//  5      OP
//  4      PLUS
//  3      COMPLEX
//  2      TIMES
//  1      COMPLEX        (pointer to A)       3 + 2i
//  0      TIMES          (pointer to B)       1 + i
//
//  (`stack.cf' and `stack.rf' are omitted to keep the example simple.)
// 
/////////////////////////////////////////////////////////////////////////////
{
  if( stack.com[0] == UNINITIALIZED )
    error("Attempt to multiply an uninitialized operator by a scalar.");
  Operator Z;
  offsetCopyStack( Z, 2, 0, 1 );
  Z.stack.com[1] = COMPLEX;
  Z.stack.com[0] = TIMES;
  Z.stack.c[0] = alpha;
  return Z;
}

Operator Operator::operator*(double a) const
//
// Scalar multiplication of an operator by a real number.
{
  if( stack.com[0] == UNINITIALIZED )
    error("Attempt to multiply an uninitialized operator by a scalar.");
  Operator Z;
  offsetCopyStack( Z, 2, 0, 0, 1 );
  Z.stack.com[1] = REAL;
  Z.stack.com[0] = TIMES;
  Z.stack.r[0] = a;
  return Z;
}

Operator Operator::operator*(ImaginaryUnit im) const
//
// Scalar multiplication of an operator by the imaginary unit `IM' or 
// by its negative `M_IM'. `IM' and `M_IM' are defined in "State.h".
{
  if( stack.com[0] == UNINITIALIZED )
    error("Attempt to multiply an uninitialized operator by a scalar.");
  Operator Z;
  offsetCopyStack( Z, 2 );
  if( im == IM )
    Z.stack.com[1] = IMAG;
  else
    Z.stack.com[1] = M_IMAG;
  Z.stack.com[0] = TIMES;
  return Z;
}

Operator Operator::operator*(ComplexFunction f) const
//
// Scalar multiplication of an operator by a complex-valued function.
{
  if( stack.com[0] == UNINITIALIZED )
    error("Attempt to multiply an uninitialized operator by a function.");
  Operator Z;
  offsetCopyStack( Z, 2, 0, 0, 0, 1 );
  Z.stack.com[1] = CFUNC;
  Z.stack.com[0] = TIMES;
  Z.stack.cf[0] = f;
  return Z;
}

Operator Operator::operator*(RealFunction f) const
//
// Scalar multiplication of an operator by a real-valued function.
{
  if( stack.com[0] == UNINITIALIZED )
    error("Attempt to multiply an uninitialized operator by a function.");
  Operator Z;
  offsetCopyStack( Z, 2, 0, 0, 0, 0, 1 );
  Z.stack.com[1] = RFUNC;
  Z.stack.com[0] = TIMES;
  Z.stack.rf[0] = f;
  return Z;
}

Operator Operator::pow(int n) const
///////////////////////////////////////////////////////////////////////
// Integer power (n > 0) of an operator.
//
// Example: If the stack of the `Operator X = A + B;' is given by
//
//  i    X.stack.com[i]    X.stack.op[i]     X.stack.c[i]     X.stack.r[i]
// 
//  2      OP
//  1      OP              (pointer to A)
//  0      PLUS            (pointer to B)
//  
// then the stack of the `Operator Z = X.pow(3);' is given by
//
//  i    Z.stack.com[i]    Z.stack.op[i]     Z.stack.c[i]     Z.stack.r[i]
//
//  10     OP
//  9      OP
//  8      PLUS
//  7      OP
//  6      OP
//  5      PLUS            (pointer to A)
//  4      TIMES           (pointer to B)
//  3      OP              (pointer to A)
//  2      OP              (pointer to B)
//  1      PLUS            (pointer to A)
//  0      TIMES           (pointer to B)
//
//  (`stack.cf' and `stack.rf' are omitted to keep the example simple.)
// 
// Each of the statements
// psi *= Z;
// or
// psi *= (A+B).pow(3);
// leads to the execution of the following code:
// {
//   State psi1(psi);
//   psi1 *= B;
//   psi *= A;
//   psi += psi1;
// }     // Here psi1 goes out of scope and is destroyed.
// {
//   State psi1(psi);
//   psi1 *= B;
//   psi *= A;
//   psi += psi1;
// }     // Here psi1 goes out of scope and is destroyed.
// {
//   State psi1(psi);
//   psi1 *= B;
//   psi *= A;
//   psi += psi1;
// }     // Here psi1 goes out of scope and is destroyed.
////////////////////////////////////////////////////////////////////
{
  if(n <= 0) error("pow(n) not defined for n <= 0.");
  if( stack.com[0] == UNINITIALIZED )
    error("Attempt to compute a power of an uninitialized operator.");
  Operator Z;
  offsetCopyStack (Z, (n-1)*(stack.comSize+1), (n-1)*stack.opSize,
    (n-1)*stack.cSize,(n-1)*stack.rSize,(n-1)*stack.cfSize,(n-1)*stack.rfSize);
  for( int i=0; i<n-1; i++ ) {
    offsetCopyStack( Z, i*(stack.comSize+1)+1, i*stack.opSize, 
	   i*stack.cSize, i*stack.rSize, i*stack.cfSize, i*stack.rfSize );
    Z.stack.com[i*(stack.comSize+1)] = TIMES;
  }
  return Z;
}

Operator Operator::hc() const
//
// Hermitian conjugate of an operator.
{
  StackPtr stackPtr = {0,0,0,0,0,0};   // Initialize stack pointers.
  return dagger( stackPtr );           // Enter the recursive `dagger'.
}

Operator Operator::dagger( StackPtr& stackPtr ) const
{
  Command c=stack.com[stackPtr.com++];   // Take command-stack content and pop.
  switch (c) 
  {
  case OPERATOR: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::dagger for case OPERATOR entered." << std::endl;
    #endif
    Operator Z;
    Z.deallocate();
    Z.allocate(1,1);
    Z.stack.com[0] = OPERATOR_HC;
    Z.stack.op[0] = stack.op[stackPtr.op++];
    return Z;
  }
  case OPERATOR_HC: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::dagger for case OPERATOR_HC entered." << std::endl;
    #endif
    Operator Z;
    Z.deallocate();
    Z.allocate(1,1);
    Z.stack.com[0] = OPERATOR;
    Z.stack.op[0] = stack.op[stackPtr.op++];
    return Z;
  }
  case COMPLEX: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::dagger for case COMPLEX entered." << std::endl;
    #endif
    Operator Z;
    Z.deallocate();
    Z.allocate(1,0,1);
    Z.stack.com[0] = COMPLEX;
    Z.stack.c[0] = conj( stack.c[stackPtr.c++] );
    return Z;   // `Z' is not a well-formed operator, but the final result is.
  }
  case REAL: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::dagger for case REAL entered." << std::endl;
    #endif
    Operator Z;
    Z.deallocate();
    Z.allocate(1,0,0,1);
    Z.stack.com[0] = REAL;
    Z.stack.r[0] = stack.r[stackPtr.r++];
    return Z;
  }
  case IMAG: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::dagger for case IMAG entered." << std::endl;
    #endif
    Operator Z;
    Z.deallocate();
    Z.allocate(1);
    Z.stack.com[0] = M_IMAG;
    return Z;
  }
  case M_IMAG: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::dagger for case M_IMAG entered." << std::endl;
    #endif
    Operator Z;
    Z.deallocate();
    Z.allocate(1);
    Z.stack.com[0] = IMAG;
    return Z;
  }
  case CFUNC: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::dagger for case CFUNC entered." << std::endl;
    #endif
    Operator Z;
    Z.deallocate();
    Z.allocate(1,0,0,0,1);
    Z.stack.com[0] = CCFUNC;
    Z.stack.cf[0] = stack.cf[stackPtr.cf++];
    return Z;   // `Z' is not a well-formed operator, but the final result is.
  }
  case CCFUNC: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::dagger for case CCFUNC entered." << std::endl;
    #endif
    Operator Z;
    Z.deallocate();
    Z.allocate(1,0,0,0,1);
    Z.stack.com[0] = CFUNC;
    Z.stack.cf[0] = stack.cf[stackPtr.cf++];
    return Z;   // `Z' is not a well-formed operator, but the final result is.
  }
  case RFUNC: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::dagger for case RFUNC entered." << std::endl;
    #endif
    Operator Z;
    Z.deallocate();
    Z.allocate(1,0,0,0,0,1);
    Z.stack.com[0] = RFUNC;
    Z.stack.rf[0] = stack.rf[stackPtr.rf++];
    return Z;   // `Z' is not a well-formed operator, but the final result is.
  }
  case PLUS: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::dagger for case PLUS entered." << std::endl;
    #endif
    Operator Z = dagger( stackPtr );
    Operator X = dagger( stackPtr );
    return X += Z;
  }
  case MINUS: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::dagger for case MINUS entered." << std::endl;
    #endif
    Operator Z = dagger( stackPtr );
    Operator X = dagger( stackPtr );
    return X -= Z;
  }
  case TIMES: {
    #ifdef DEBUG_TRACE
    std::cout << "Operator::dagger for case TIMES entered." << std::endl;
    #endif
    Operator Z = dagger( stackPtr );
    Operator X = dagger( stackPtr );
    return Z * X;     // Notice the reversal of the order of terms.
  }
  case UNINITIALIZED:    
    error("Attempt to compute the conjugate of an uninitialized operator.");

  default:
    error("dagger: Unknown object in command stack.");
  }
}

//////////////////////////////////////////////////////////////////////
// The following algebraic operations on operators are all defined in
// terms of algebraic operations defined above.
//////////////////////////////////////////////////////////////////////

Operator operator-(const Operator& op)      // Unary `-'.
{
  Operator result = -1.0 * op;
  return result;
}

Operator operator+(const Operator& op)      // Unary `+'.
{
  Operator result = op;
  return result;
}

Operator& Operator::operator+=(const Operator& op)
{
  *this = *this + op;
  return *this;
}

Operator& Operator::operator-=(const Operator& op)
{
  *this = *this - op;
  return *this;
}

Operator& Operator::operator*=(const Complex& alpha)
{
  *this = *this * alpha;
  return *this;
}

Operator& Operator::operator*=(double a)
{
  *this = *this * a;
  return *this;
}

Operator& Operator::operator*=(ImaginaryUnit im)
{
  *this = *this * im;
  return *this;
}

Operator operator*(const Complex& alpha, const Operator& op)
{
  return op * alpha;
}
  
Operator operator*(double a, const Operator& op)
{
  return op * a;
}
  
Operator operator*(ImaginaryUnit im, const Operator& op)
{
  return op * im;
}

Operator& Operator::operator*=(ComplexFunction f)
{
  *this = *this * f;
  return *this;
}

Operator& Operator::operator*=(RealFunction f)
{
  *this = *this * f;
  return *this;
}

Operator operator*(ComplexFunction f, const Operator& op)
{
  return op * f;
}
  
Operator operator*(RealFunction f, const Operator& op)
{
  return op * f;
}
