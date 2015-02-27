//   State.h -*- C++ -*- State algebra in Hilbert space.
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

#ifndef _State_hhh
#define _State_hhh 1

#include "Complex.h"

enum FreedomType{ ALL, FIELD, SPIN, ATOM };
// FreedomType lists the different types of physical systems specifically
// recognized by this program; more may be added if necessary.
// FIELD is an oscillator degree of freedom, which is described in the
//   Fock state or Excited Coherent State basis;
// SPIN is a spin 1/2 or two-level atom;
// ATOM is an atom or N-level system;
// ALL leaves the physical type of system unspecified.
// Note that certain member functions will only work for degrees of freedom
// of a particular type.

enum ImaginaryUnit{ IMAGINARY_UNIT, MINUS_IMAGINARY_UNIT };
const ImaginaryUnit IM = IMAGINARY_UNIT;
const ImaginaryUnit M_IM = MINUS_IMAGINARY_UNIT;
// For efficiency, special routines have been included allowing States
// and Operators to be multiplied by i and -i without invoking the full
// Complex arithmetic.  Note that multiplication by ImaginaryUnit is
// define ONLY for States and Operators; 3*IM is not a legal
// expression.

#ifndef NON_GNU_PROTOTYPE
  class PrimaryOperator;
#else
  extern class PrimaryOperator;
#endif

// I/O functions for FreedomType

ostream& operator<<( ostream&, FreedomType );
istream& operator>>( istream&, FreedomType& );

class State{
// The State class represents quantum states in a particular choice of
// Hilbert space; this can include varying numbers of physical degrees
// of freedom, described by varying numbers of basis states.

public:						// public functions

// constructors and destructors

  State();
    // Default constructor; produces a State of size 0.
  State(const State& state);
    // Copy constructor.
  State(int n, FreedomType form=FIELD);
    // Produces a single degree-of-freedom State with n basis states;
    // all basis states have amplitude 0 except the ground state,
    // which has amplitude 1; form gives the FreedomType of the State
    // (default is FIELD).
  State(int n, int* dimensions, FreedomType* forms);
    // Multiple degree of freedom ground State.  n gives the number of
    // degrees of freedom; dimensions is an array of n integers, which
    // specify the number of basis states to allocate for each degree
    // of freedom; forms is an array of n FreedomTypes, giving the
    // physical type of each degree of freedom.  The basis of the full
    // n-freedom state is given by the products of the basis states of
    // the n freedoms; all have amplitude 0 except for the ground state.
  State(int n, Complex* elements, FreedomType form=FIELD);
    // Produces a one degree-of-freedom State with n basis
    // states.  elements is an array of n Complex numbers, representing
    // the amplitudes of the n basis states; form gives the FreedomType
    // of the State (default is FIELD).
  State(int n, int nstate, FreedomType form=FIELD);	// Fock state
    // Produces a single degree-of-freedom State with n basis states;
    // all basis states have amplitude 0 except state number nstate,
    // which has amplitude 1; form gives the FreedomType of the State
    // (default is FIELD).
  State(int n, Complex alpha, FreedomType form=FIELD);	// Coherent state
    // Produces a one degree-of-freedom State with n basis
    // states in a coherent state given by the Complex number alpha.
    // The State is represented in a Fock (number) state basis.
    // form gives the FreedomType of the State (default is FIELD).
  State(int n, int nstate, Complex alpha, FreedomType form=FIELD);
    // Excited coherent state.  Produces a single degree-of-freedom State
    // with n basis states.  Basis state number nstate has amplitude 1,
    // all others have amplitude 0, and the state is in the excited coherent
    // state or displaced Fock state basis, centered at alpha in phase space.
    // form gives the FreedomType of the State (default is FIELD).
  State(int n, State* stateList);		// Product state
    // Produces an n degree-of-freedom state.  stateList is an array of
    // n one-freedom states.  The n-freedom state will be produced in a
    // product state of the n states in stateList.  The FreedomTypes and
    // dimensions of the different freedoms of the n-freedom state will
    // match those of the one-freedom states in stateList.
  ~State();					// destructor

// public functions used by constructors and destructors

  void fock(int n, int nstate);			// create a Fock state
  void coherent(int n, Complex alpha);		// create a Coherent state
  void productState(int n, State* stateList);	// create a product state

// Member arithmetic operations

  inline Complex& operator[](int n) {	        // subscript operator; gives
     if( nSkip != 1 )				// access to the amplitude
       return myPointer[nSkip*n];		// of the nth basis state
     else					// Note that this amplitude
       return myPointer[n];			// can be changed as well as
  };						// read.  Usage: psi[n]
						// Note the presence of nSkip;
						// this is part of the
						// SkipVector structure, used
						// for multiple freedom states
  Complex elem(const int*) const;
    // MultiDim subscripting; takes as an argument an array of integers
    // of length equal to the number of degrees of freedom.  Returns
    // the amplitude of the corresponding basis state.  Note that this
    // subscripting is read-only.  Usage: psi.elem(n_array)
  Complex& operator[](int*);
    // same as elem (above), but also permits the amplitudes to be
    // changed.  Usage:  psi[n_array]
  State& operator=(const State&);		// assignment
  State& operator=(int);
    // zero assignment; enables one to type psi=0 to set all amplitudes
    // to 0.  Gives an error for any int other than 0.
  Complex operator*(const State&) const;	// inner product
  State& operator*=(const Complex&);		// multiply by Complex scalar
  State& operator*=(double);			// multiply by real scalar
  State& operator*=(ImaginaryUnit);		// multiply by ImaginaryUnit
  State& operator+=(const State&);		// add a State
  State& operator-=(const State&);		// subtract State

// Friend arithmetic operations

  friend State operator*(const Complex&, const State&);
    // multiply a State by a Complex scalar:  z*psi
  friend State operator*(const State&, const Complex&);
    // multiply a State by a Complex scalar (other order): psi*z
  friend State operator*(double, const State&);
    // multiply a State by a real scalar:  x*psi
  friend State operator*(const State&, double);
    // multiply a State by a real scalar (other order):  psi*x
  friend State operator*(ImaginaryUnit, const State&);
    // multiply a State by an ImaginaryUnit (i or -i)
  friend State operator*(const State&, ImaginaryUnit);
    // multiply a State by an ImaginaryUnit (i or -i) (other order)
  friend State operator+(const State&, const State&);
    // add two States
  friend State operator-(const State&, const State&);
    // subtract one State from another
  friend State operator+(const State&);		// unary +
  friend State operator-(const State&);		// unary -

// Friend I/O operations

  friend ostream& operator<<( ostream&, const State&);
    // outputs a state in a standard ASCII form.  This can be used to
    // save and recover results of a calculation.
  friend istream& operator>>( istream&, State& );
    // inputs a state in a standard ASCII form.  This can be used to
    // save and recover results of a calculation.

// Information-returning and utility member functions

  void xerox(const State& a);
    // make MINIMAL copy of State (for use in temps).  Improves efficiency
    // when dynamical allocation of basis states is being used.  Chiefly
    // used by the Operator class; should not be needed by ordinary users.
    // For an explanation of the dynamical allocation see adjustCutoff
    // below.
  int size();			// length of data array
  int getSize(int = 0);		// size of nth degree of freedom
  void diagnostic();		// debugging info
  void normalize();		// normalize state, i.e., psi*psi = 1

// Member functions accessing coordinates; for a full explanation of
// this, see the basis-changing member functions below.

  Complex centerCoords();
    // return center of coordinates
  Complex getCoords(int = 0);
    // center of coordinates of nth freedom (default 0)
  void setCoords(Complex&,int=0);
    // set the value of the coords for a freedom (default 0)
  void displaceCoords(Complex&, int=0);
    // adds a Complex displacement to the center of
    // coordinates of a freedom (default freedom is 0)
  double checkBounds(int, int=2);
    // check amplitudes of top basis states
    // (default: top 2 basis states)

// Basis-changing member functions

// This QSD library makes use of the localization property to greatly
// improve the efficiency of calculations.  In QSD, for a wide variety
// of problems, field states tend towards highly localized wavepackets
// in phase space.  These wave packets are closely centered on some point
// alpha (a Complex number) which is given by the expectation value
// alpha = <a>, where a is the harmonic oscillator annihilation operator.
//
// If alpha is large, it requires a great many ordinary Fock states to
// represent such a wavepacket; the number of Fock states n goes like
// |alpha|^2.  By choosing a different set of basis states an enormous
// savings is possible.
//
// We use the excited coherent state basis |alpha,n> to represent our
// states, choosing alpha=<a> for maximum efficiency.  (Ordinary Fock
// states would correspond to alpha=0.)  As the value of <a> will change
// with time, the basis must also be changed fairly often.  This adds to
// the cost of a calculation; but the savings from the moving basis
// far outweigh this added complexity.
//
// Note that this only applies to freedoms of type FIELD.  For multiple
// degree-of-freedom states, each FIELD degree of freedom can be moved
// separately.  Trying to move the basis of a non-FIELD freedom will
// produce an error.
//
// Note also that the user is not required to use the moving basis.  For
// problems without strong localization, the moving basis adds to the
// computational overhead while producing little benefit, and should not
// be invoked.
//
// For further details, see J. Phys. A 28, 5401-5413 (1995).

  void moveCoords(const Complex& displacement, int theFreedom=0,
      double shiftAccuracy=1e-4);
    // Relative shift of the center of coordinates.  displacement
    // gives the amount by which to shift alpha, theFreedom indicates
    // which degree of freedom is to be shifted.  shiftAccuracy
    // gives the accuracy with which to make the shift
    // (default 1e-4).  The physical state is unchanged, but it is
    // represented in a new basis |alpha+displacement,n>
    // Uses private moveStep member function.
  void recenter(int theFreedom=0,double shiftAccuracy=1e-4);
    // Recenters freedom theFreedom at its expectation value in phase space.
    // The physical state is unchanged, but is represented in a new
    // basis |<a>,n>, a being the annihilation operator for the
    // selected degree of freedom.  Uses moveCoords.
  void moveToCoords(const Complex& alpha, int theFreedom=0,
      double shiftAccuracy=1e-4);
    // Like moveCoords, but instead of shifting by a Complex displacement
    // it moves the basis to a new absolute position alpha in phase space.
    // Uses moveCoords.
  void centerOn(State& psi,double shiftAccuracy=1e-4);
    // Leaves the physical state unchanged, but represents it in the
    // same basis as the given state psi.  Uses moveCoords.  Can be used
    // with multiple degree-of-freedom states, though the states must
    // have the same number and type of freedoms; it will change the
    // basis only of the FIELD degrees of freedom.

// Cutoff-adjusting member functions

// For efficiency, it is possible to restrict the number of basis vectors
// actually used by including only those with amplitudes appreciably
// greater than zero.  The criterion used is based on two parameters:
// epsilon and padSize.  All of the top states of a freedom whose total
// probability is less than epsilon are excluded except for a number of
// "buffer" basis states, or "pad", equal to padSize.  If the "pad" begins
// to gain appreciable probability (i.e., probability over epsilon) the
// number of basis states can be dynamically increased.  In this way,
// no more basis states are used than are necessary.  This is particularly
// useful in conjunction with the moving basis.

  void adjustCutoff(int theFreedom=0, double epsilon=1e-4, int padSize=2);
    // adjust amount of storage used by the State by adjusting the
    // cutoff for freedom number theFreedom.
  void fullSize();	// makes size of state match physical size of storage

private:			// private data and functionss

// SkipVector part -- used to act on single degree of freedom in memory
//
// The QSD code assumes that all Operators are defined in terms of Primary
// Operators which act on a single degree of freedom.  In order for this
// to work successfully, it must be possible to loop over a single degree
// of freedom, leaving all the others unchanged.  This is embodied in the
// notion of a SkipVector -- a data structure which is treated like an
// ordinary array, but steps through memory by an ``skip'' larger than one.

  Complex* myPointer;		// pointer to storage
  int mySize;			// array size
  int nSkip;			// steps between elements (the ``skip'')

// One dimensional state part --
// used for greater efficiency in 1 Freedom problems

  int totalDim;			// total number of elements
  int maxSize;			// number of elements actually used
  FreedomType myType;		// type of freedom
  Complex* data;		// element storage; this points to an array
				// of complex numbers giving the amplitudes
				// of the basis states
  Complex coord;		// center of coordinates in phase space
				// (see the moving basis, above)

// MultiState part -- used for >1 degrees of freedom
//
// The basis states for a multiple-degree-of-freedom state are the
// products of the basis states of the individual degrees of freedom
// which make up the total state.  If there are n_i basis states to
// represent the ith degree of freedom, then a particular basis state
// of the total state is given by indices {i_0,...,i_N} for an N+1 freedom
// state.
//
// The amplitudes for these basis states are stored in the array data,
// just as for single degree-of-freedom states.  The amplitude for the
// basis state {i_0,...,i_N} is stored at the location
//
//   loc = i_0 + i_1*n_0 + i_2*n_0*n_1 + ... + i_N*n_0*...*n_(N-1).
//
// where i_m ranges from 0 to n_m-1.
//
// Incrementing the mth index by 1 is equivalent to stepping through
// the data array by an amount
//
//   nSkips[m] = n_0*n_1*...*n_(m-1).

  int nFreedoms;		// number of degrees of freedom
  int* nDims;			// number of dimensions (basis states)
				// for each freedom
  int* sizes;			// number of dimensions (basis states)
				// actually used by each freedom
  int* nSkips;			// index spacing in data field of ith freedom
  int* partDims;		// subspace dim of first i freedoms
  FreedomType* freedomTypes;	// types of degrees of freedom (FIELD, SPIN, etc.)
  Complex* coords;		// centers of coordinates in phase space

// Private member functions

  void apply(PrimaryOperator&, int, int, FreedomType, double);
    // apply a PrimaryOperator to a State
  void copy(const State&);		// copy a state
  void free();				// recycle storage
  void moveStep(Complex&, Complex&);	// shift coords by infinitesimal step
  void stretchFreedom(int,int);		// increase storage used by a freedom
  void shrinkFreedom(int,int);		// compact storage used by a freedom
  void error(const char*) const;	// print out error message and exit

  friend class Operator;			// friend class
};

#endif
