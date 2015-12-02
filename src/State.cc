//   State.cc -- State algebra in Hilbert space.
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
//   Dr. Todd Brun              Tel    +44 (0)171 775 3292
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
#include <math.h>
#include <cstdlib>

#include "State.h"
#include "PrimOp.h"

static const char rcsid[] = "$Id: State.cc,v 3.1 1996/11/19 10:05:08 rschack Exp $";

// These are global variables of file scope, used to avoid recalculating
// the square roots of the integers more often than necessary.

static double smallSize = 1.0e-8;
static int intSqrtMax = 0;
static double* intSqrtList = 0;

// This utility function is used to maintain the list of square roots.
// This function has file scope.

static void setIntSqrtMax(int maximum)
{
  int newsize,i;
  double* newPtr;

  if( maximum > intSqrtMax) {      // make sure of size
    newsize = 2*maximum;      // allow plenty of room
    newPtr = new double[newsize];    // allocate space
    for (i=0; i<intSqrtMax; i++)
      newPtr[i] = intSqrtList[i];    // re-use calculated values
    for (i=intSqrtMax; i<newsize; i++)
      newPtr[i] = sqrt(i);      // calculate new values
#ifndef NON_GNU_DELETE
    if(intSqrtList != 0) delete [] intSqrtList;  // recycle memory
#else // ---- NON GNU CODE ----
    if(intSqrtList != 0) delete [intSqrtMax] intSqrtList; // recycle memory
#endif // ---- NON GNU CODE ----
    intSqrtList = newPtr;      // point to new list
    intSqrtMax = newsize;      // remember new size
  }
}

//  I/O operators for enum FreedomType; necessary of I/O of States

std::ostream& operator<<( std::ostream& s, FreedomType theType )
  //
  // Output a FreedomType to an output stream
{
  switch(theType) {
    case ALL:
      s << 0;
      break;
    case FIELD:
      s << 1;
      break;
    case SPIN:
      s << 2;
      break;
    case ATOM:
      s << 3;
      break;
    default:
      std::cerr << "Unknown FreedomType outputted in <<!" << std::endl;
      exit(1);
  }
  return s;
}

std::istream& operator>>( std::istream& s, FreedomType& theType )
  //
  // Input a FreedomType from an input stream
{
  int i;

  s >> i;
  switch(i) {
    case 0:
      theType = ALL;
      break;
    case 1:
      theType = FIELD;
      break;
    case 2:
      theType = SPIN;
      break;
    case 3:
      theType = ATOM;
      break;
    default:
      std::cerr << "Unknown FreedomType inputted in >>!" << std::endl;
      exit(1);
  }
  return s;
}

// copy and free utility functions used by Constructors, Destructors, and
// assignment operator.

void State::free() {
  if( data != 0) {      // if memory is allocated...
#ifndef NON_GNU_DELETE
    delete [] data;      // recycle storage
    if(nFreedoms > 1) {
      delete [] nDims;      // recycle parameter arrays
      delete [] sizes;
      delete [] nSkips;
      delete [] partDims;
      delete [] freedomTypes;
      delete [] coords;
    }
#else // ---- NON GNU CODE ----
    delete [totalDim] data;    // NON GNU recycle storage
    if( nFreedoms > 1 ) {
      delete [nFreedoms] nDims;    // NON GNU recycle parameter arrays
      delete [nFreedoms] sizes;
      delete [nFreedoms] nSkips;
      delete [nFreedoms] partDims;
      delete [nFreedoms] freedomTypes;
      delete [nFreedoms] coords;
    }
#endif // ---- NON GNU CODE ----
    nFreedoms = 0;      // zero variables
    totalDim = 0;
    maxSize = 0;
    myType = ALL;
    mySize = 0;
    nSkip = 0;
    coord = 0;
    myPointer = 0;      // reset pointers to zero
    data = 0;
    nDims = 0;
    sizes = 0;
    nSkips = 0;
    coords = 0;
    partDims = 0;
    freedomTypes = 0;
  }
}

void State::copy(const State& a)  // make a copy of an existing state
{
  int i;
  if( nFreedoms != a.nFreedoms) {  // if no memory allocated, allocate it
    nFreedoms = a.nFreedoms;
    if( nFreedoms > 1 ) {
      nDims = new int[nFreedoms];
      sizes = new int[nFreedoms];
      nSkips = new int[nFreedoms];
      partDims = new int[nFreedoms];
      freedomTypes = new FreedomType[nFreedoms];
      coords = new Complex[nFreedoms];
    }
  }
  if( totalDim != a.totalDim ) {  // if no memory allocated, allocate it
    totalDim = a.totalDim;
    data = new Complex[totalDim];
  }
  maxSize = a.maxSize;      // copy member data
  nSkip = a.nSkip;
  coord = a.coord;
  mySize = a.mySize;
  myType = a.myType;
  myPointer = data;
  for (i=0; i<maxSize; i++)
    data[i] = a.data[i];
  if( nFreedoms > 1 ) {
    for (i=0; i<nFreedoms; i++) {
      nDims[i] = a.nDims[i];
      sizes[i] = a.sizes[i];
      nSkips[i] = a.nSkips[i];
      partDims[i] = a.partDims[i];
      freedomTypes[i] = a.freedomTypes[i];
      coords[i] = a.coords[i];
    }
  }
}

void State::xerox(const State& a)
  // Unlike copy, xerox makes a copy of an existing state including only
  // the memory which is actually being used, rather than all the memory
  // allocated.  This is useful as a more efficient way of creating
  // temporary states than using the normal copy-constructor.
  //
  // Note that xerox should NOT be used to assign to a state which already
  // has memory allocated.  This results in inefficient memory handling, and
  // might cause errors.
{
  int i;
#ifndef OPTIMIZE_QSD
  if( data != 0 )
    error("State::xerox should not be used on a state with allocated memory!");
#endif
  nFreedoms = a.nFreedoms;
  if( nFreedoms > 1 ) {
    nDims = new int[nFreedoms];
    sizes = new int[nFreedoms];
    nSkips = new int[nFreedoms];
    partDims = new int[nFreedoms];
    freedomTypes = new FreedomType[nFreedoms];
    coords = new Complex[nFreedoms];
  }
  totalDim = a.maxSize;
  data = new Complex[totalDim];
  maxSize = a.maxSize;      // copy member data
  nSkip = a.nSkip;
  coord = a.coord;
  mySize = a.mySize;
  myType = a.myType;
  myPointer = data;
  for (i=0; i<maxSize; i++)
    data[i] = a.data[i];
  if( nFreedoms > 1 ) {
    for (i=0; i<nFreedoms; i++) {
      nDims[i] = a.sizes[i];
      sizes[i] = a.sizes[i];
      nSkips[i] = a.nSkips[i];
      partDims[i] = a.partDims[i];
      freedomTypes[i] = a.freedomTypes[i];
      coords[i] = a.coords[i];
    }
  }
}

State::State() {  // Constructor for uninitialized State
  nFreedoms = 0;
  totalDim = 0;
  maxSize = 0;
  nDims = 0;
  sizes = 0;
  nSkips = 0;
  partDims = 0;
  freedomTypes = 0;
  coords = 0;
  data = 0;
  myPointer = 0;
  nSkip = 0;
  myType = ALL;
  mySize = 0;
  coord = 0;
}

State::State(int n, FreedomType form)     // Constructor for 1D
{            // ground state
#ifndef OPTIMIZE_QSD
  if (n<2) error("Invalid dimension size in State constructor!");
  if( (form==SPIN) && (n != 2) ) error("Spins must have dimension 2!");
#endif
  nFreedoms = 1;
  totalDim = maxSize = mySize = n;
  nSkip = 1;
  coord = 0;
  myType = form;
  data = new Complex[totalDim];
  myPointer = data;
  nDims = 0;
  sizes = 0;
  nSkips = 0;
  partDims = 0;
  freedomTypes = 0;
  coords = 0;
  data[0] = 1.0;        // put in the ground state
  for (int i=1; i<totalDim; i++) data[i] = 0;
}

State::State(int n, int* dimensions, FreedomType* forms) // Constructor for
{               // multidimensional
  int i;             // ground state
#ifndef OPTIMIZE_QSD
  if (n<1) error("Illegal number of freedoms in State constructor!");
  if (n==1) error("Do not call this constructor for 1-freedom states!");
#endif
  nFreedoms = n;
  nDims = new int[n];      // allocate memory
  sizes = new int[n];
  nSkips = new int[n];
  partDims = new int[n];
  freedomTypes = new FreedomType[n];
  coords = new Complex[n];
  totalDim = 1;
  for (i=0; i<n; i++) {
#ifndef OPTIMIZE_QSD
    if(dimensions[i] < 1) error("Illegal dimension in State constructor!");
    if( (forms[i] == SPIN) && (dimensions[i] != 2) )
      error("Spin states must have 2 dimensions!");
#endif
    nSkips[i] = totalDim;
    totalDim *= (nDims[i] = dimensions[i]);
    sizes[i] = nDims[i];
    partDims[i] = totalDim;
    freedomTypes[i] = forms[i];
    coords[i] = 0;
  }
  maxSize = totalDim;
  data = new Complex[totalDim];    // allocate memory
  myPointer = data;
  mySize = totalDim;
  nSkip = 1;
  coord = 0;
  myType = freedomTypes[0];
  data[0] = 1;
  for (i=1; i<totalDim; i++) data[i] = 0;
}

State::State(int n, int nstate, FreedomType form)  // Constructor for
{              // Fock state
#ifndef OPTIMIZE_QSD
  if (n<2) error("Invalid dimension size in State constructor!");
  if( form != FIELD )
    error("Constructor only valid for FIELD type.");
#endif
  data = 0;
  fock(n,nstate);    // call Fock state function
}

State::State(int n, Complex alpha, FreedomType form)     // Constructor for
{               // coherent state
#ifndef OPTIMIZE_QSD
  if (n<2) error("Invalid dimension size in State constructor!");
  if( form != FIELD )
    error("Constructor only valid for FIELD type.");
#endif
  data = 0;
  coherent(n,alpha);    // call coherent state function
}

State::State(int n, int nstate, Complex alpha, FreedomType form)
  //
  // Constructor for excited coherent state
{
#ifndef OPTIMIZE_QSD
  if (n<2) error("Invalid dimension size in State constructor!");
  if( form != FIELD )
    error("Constructor only valid for FIELD type.");
#endif
  data = 0;
  fock(n,nstate);  // call Fock state function
  coord = alpha;  // displace in phase space
}

State::State(int n, Complex* elements, FreedomType form)
  //
  // Constructor for general 1 degree of freedom state
{
#ifndef OPTIMIZE_QSD
  if (n<2) error("Invalid dimension size in State constuctor!");
  if( (form==SPIN) && (n != 2) ) error("Spins must have dimension 2!");
#endif
  nFreedoms = 1;
  totalDim = maxSize = mySize = n;
  nSkip = 1;
  coord = 0;
  myType = form;
  data = new Complex[totalDim];
  myPointer = data;
  nDims = 0;
  sizes = 0;
  nSkips = 0;
  partDims = 0;
  freedomTypes = 0;
  coords = 0;
  for (int i=0; i<totalDim; i++) data[i] = elements[i];
}

State::State(const State& a)       // Copy-initializer
{
  nFreedoms = 0;
  totalDim = 0;
  copy(a);
}

State::State(int n, State* stateList)    // Constructor for
{            // product state
#ifndef OPTIMIZE_QSD
  if (n<1) error("Illegal number of freedoms in State constructor!");
#endif
  data = 0;
  productState(n,stateList);   // call product state function
}

State::~State() { free(); }    // Destructor

void State::fock(int n, int nstate)  // Create a Fock State
{
#ifndef OPTIMIZE_QSD
  if( (nstate >= n) || (n < 2) || (nstate < 0) )  // check dimension
    error("Illegal arguments to State::fock!");    // exit on error
#endif
  free();        // recycle memory
  totalDim = maxSize = mySize = n;  // set up state
  nFreedoms = 1;
  nSkip = 1;
  coord = 0;
  myType = FIELD;
  data = new Complex[totalDim];    // allocate memory
  myPointer = data;
  nDims = 0;
  sizes = 0;
  nSkips = 0;
  partDims = 0;
  freedomTypes = 0;
  coords = 0;
  for (int i=0; i<totalDim; i++) data[i] = 0;
  data[nstate] = 1;          // Fock state
}

void State::coherent(int n, Complex alpha)  // Create a coherent state
{
#ifndef OPTIMIZE_QSD
  if( n < 2 )
    error("Invalid state size in State::coherent!");   // exit on error
  if( norm(alpha) > n )
    error("Coherent state too big for storage in State::coherent!");
#endif
  Complex z;
  free();        // recycle memory
  totalDim = maxSize = mySize = n;  // set up state
  nFreedoms = 1;
  nSkip = 1;
  coord = 0;
  myType = FIELD;
  data = new Complex[totalDim];    // allocate memory
  myPointer = data;
  nDims = 0;
  sizes = 0;
  nSkips = 0;
  partDims = 0;
  freedomTypes = 0;
  coords = 0;
  z = exp( - 0.5*norm(alpha));    // normalization factor
  for (int i=0; i<totalDim; i++) {
    data[i] = z;      // coherent state elements
    z *= alpha/sqrt(i+1);
  }
}

void State::productState(int n, State* stateList)  // Create a product state
{
  int i,j,k,l;
  Complex element;
#ifndef OPTIMIZE_QSD
  if( n < 1 ) error("Invalid number of degrees of freedom in productState!");
#endif
  free();
  nFreedoms = n;
  nDims = new int[n];      // allocate memory
  sizes = new int[n];
  nSkips = new int[n];
  partDims = new int[n];
  freedomTypes = new FreedomType[n];
  coords = new Complex[n];
  totalDim = 1;
  for (i=0; i<n; i++) {      // set up memory structure for data
#ifndef OPTIMIZE_QSD
    if( stateList[i].nFreedoms != 1 )
      error("Can only produce products of 1D states!");
    if( stateList[i].totalDim < 2 )
      error("Invalid dimension size in productState!");
#endif
    nSkips[i] = totalDim;
    totalDim *= (nDims[i] = stateList[i].totalDim);
    sizes[i] = nDims[i];
    partDims[i] = totalDim;
    freedomTypes[i] = stateList[i].myType;
    coords[i] = stateList[i].coord;
  }
  mySize = totalDim;
  maxSize = totalDim;
  nSkip = 1;        // initialize state variables
  coord = coords[0];
  myType = freedomTypes[0];
  data = new Complex[totalDim];
  myPointer = data;
  for (i=0; i<totalDim; i++) {
    k = i;
    element = 1;
    for (j=n-1; j>=0; j--) {
      l = k/nSkips[j];      // index i_j of jth state
      element *= stateList[j][l];  // product state
      if( element == 0 ) break;    // if zero, leave this inner loop
      k %= nSkips[j];
    }
    data[i] = element;  // psi_1[i_1]*psi_2[i_2]*psi_3[i_3]*...*psi_n[i_n]
  }
}

State& State::operator=(const State& a)    // assignment operator
{
  if (this != &a) {  // otherwise a=a would lead to trouble
    if( (totalDim != a.totalDim) || (nFreedoms != a.nFreedoms) )
      free();    // recycle any memory
    copy(a);    // copy a into b
  }
  return *this;
}

State& State::operator=(int n)    // assign zero operator
{
  if (n != 0) error("Cannot assign a State a value other than zero!");
  for (int i=0; i<maxSize; i++) data[i] = 0;
  return *this;
}

Complex& State::operator[](int* indexList)
  //
  // Multiple degree of freedom subscript operator
  // DO NOT USE THIS WITH UNNAMED TEMPORARIES
  //
{
  int i,j=0;
  if (nFreedoms == 1)
    return data[indexList[0]];    // 1 degree of freedom case...
  for (i=0; i<nFreedoms; i++) j += indexList[i]*nSkips[i];
#ifndef OPTIMIZE_QSD
  if( (j>maxSize) || (j < 0) )
    error("Illegal arguments passed to function State::operator[]!");
#endif
  return data[j];
}

Complex State::elem(const int* indexList) const
//
// Multiple degree of freedom subscript operator.  Use this function
// when writing to the state might cause an error (e.g., when reading
// an element from an unnamed temporary).
//
{
  int i,j=0;
  if (nFreedoms == 1)
    return data[indexList[0]];    // 1 degree of freedom case...
  for (i=0; i<nFreedoms; i++) j += indexList[i]*nSkips[i];
#ifndef OPTIMIZE_QSD
  if( (j>maxSize) || (j < 0) )
    error("Illegal arguments passed to function elem!");
#endif
  return data[j];
}

Complex State::operator*(const State& a) const    // inner product
{
#ifndef OPTIMIZE_QSD
  if( mySize != a.mySize )
    error("Incompatible state sizes in State::operator*!");
  if( myType != a.myType )
    error("Incompatible state types in State::operator*!");
#endif
  int i,j,k;
  Complex sum = 0.0;
  for (i=0,j=0,k=0; i<mySize; i++,j+=nSkip,k+=a.nSkip)
    sum += conj(myPointer[j])*a.myPointer[k];
  return sum;
}

State& State::operator*=(const Complex& z)  // multiplication by Complex
{
  int i,j;
  for (i=0,j=0; i<mySize; i++,j+=nSkip)
    data[j] *= z;
  return *this;
}

State& State::operator*=(double x)    // multiplication by real
{
  int i,j;
  for (i=0,j=0; i<mySize; i++,j+=nSkip)
    data[j] *= x;
  return *this;
}

State& State::operator*=(ImaginaryUnit im)  // multiplication by i or -i
{
  int i,j;
  if( im == IM )
    for (i=0,j=0; i<mySize; i++,j+=nSkip) data[j].timesI();
  else
    for (i=0,j=0; i<mySize; i++,j+=nSkip) data[j].timesMinusI();
  return *this;
}

State& State::operator+=(const State& a)  // add state
{
#ifndef OPTIMIZE_QSD
  if( mySize != a.mySize )
    error("Incompatible state sizes in State::operator+=!");
  if( myType != a.myType )
    error("Incompatible state types in State::operator+=!");
#endif
  int i,j,k;
  for (i=0,j=0,k=0; i<mySize; i++,j+=nSkip,k+=a.nSkip)
    data[j] += a.data[k];
  return *this;
}

State& State::operator-=(const State& a)  // subtract state
{
#ifndef OPTIMIZE_QSD
  if( mySize != a.mySize )
    error("Incompatible state sizes in State::operator+=!");
  if( myType != a.myType )
    error("Incompatible state types in State::operator+=!");
#endif
  int i,j,k;
  for (i=0,j=0,k=0; i<mySize; i++,j+=nSkip,k+=a.nSkip)
    data[j] -= a.data[k];
  return *this;
}

void State::normalize()      // normalize state to 1
{
  int i;
  double sum=0.0;
  double epsLimit;
  epsLimit = smallSize/maxSize;
  for (i=0; i<maxSize; i++) {
    sum += norm(data[i]);
  }
  if(sum > epsLimit ) {      // check if zero state...
    sum = 1.0/sqrt(sum);
    for (i=0; i<maxSize; i++)
      data[i] *= sum;
  }
}

double State::checkBounds(int theFreedom, int numChecks)
  //
  // Check how much probability is in the top numChecks levels
  // of the degree of freedom theFreedom.
  //
{
#ifndef OPTIMIZE_QSD
  if( (theFreedom < 0) || (theFreedom >= nFreedoms) )
    error("Illegal freedom passed to State::checkBounds!");
  if( (nFreedoms == 1) && (numChecks > maxSize) )
    error("Too many checks requested in State::checkBounds!");
  if( (nFreedoms > 1) && (numChecks > sizes[theFreedom]) )
    error("Too many checks requested in State::checkBounds!");
#endif
  double sum = 0.0;
  int i,j,k,l;
  if( nFreedoms == 1) {        // 1 freedom case
    for(i=1; i<=numChecks; i++) {
      sum += norm(data[maxSize - i]);    // add up probability
    }
  }
  else {          // multi-freedom case
    int theSkip = nSkips[theFreedom];
    int theBound = (sizes[theFreedom] - 1)*theSkip;
    int theDim = partDims[theFreedom];
    for(i=0; i<theSkip; i++)      // loop over other indices
      for(j=0; j<maxSize; j+=theDim)
        for(k=0,l=theBound; k<numChecks; k++,l-=theSkip) {
          sum += norm(data[i+j+l]);      // add up probability
        }
  }
  return sum;
}

void State::fullSize()
  //
  // Makes the cutoffs of the state match the physical size in memory
  //
{
  if( nFreedoms == 1 ) {
    for (int i=maxSize; i<totalDim; i++)
      data[i] = 0.0;
    maxSize = totalDim;
    mySize = totalDim;
  }
  else {
    for(int i=0; i<nFreedoms; i++)
      if( sizes[i] < nDims[i] ) stretchFreedom(i,nDims[i]);
  }
}

void State::adjustCutoff(int theFreedom, double epsilon, int padSize)
  //
  // This member function adjusts the amount of memory actually being
  // used by a State.  As very frequently only the lowest few states have
  // significant probability, it is more efficient not to loop over the
  // top states at all.  This routine eliminates the top N states (which
  // contain total probability less than epsilon) with a pad of size padSize
  // for safety.
  //
{
#ifndef OPTIMIZE_QSD
  if( (theFreedom < 0) || (theFreedom >= nFreedoms) )
    error("Illegal freedom passed to State::checkBounds!");
  if( (nFreedoms == 1) && (padSize > totalDim) )
    error("Too large a padSize requested in State::adjustCutoff!");
  if( (nFreedoms > 1) && (padSize > nDims[theFreedom]) )
    error("Too large a padSize requested in State::adjustCutoff!");
  if( padSize < 1 ) error("padSize must be positive in State::adjustCutoff!");
  if( (epsilon < 0.0) || (epsilon > 1.0) )
    error("Illegal value of epsilon in State::adjustCutoff!");
  if( (nFreedoms == 1) && (myType != FIELD) )
    error("Can only adjust cutoff size for fields!");
  if( (nFreedoms > 1) && (freedomTypes[theFreedom] != FIELD) )
    error("Can only adjust cutoff size for fields!");
#endif
  double sum = 0.0;
  int i,j,l;
  normalize();        // make sure total prob. is 1
  if( nFreedoms == 1) {      // 1 freedom case
    for(i=maxSize; sum < epsilon;)
      sum += norm(data[--i]);    // add up probability
    i += (padSize + 1);
    if( i > totalDim ) {    // enough states allocated?
      std::cerr << "Warning!  Significant probability in top states!" << std::endl;
      i = totalDim;
    }
    if( i > maxSize )        // stretch # of states;
      for(j=maxSize; j<i; j++) data[j] = 0;  // zero top states
    maxSize = i;
    mySize = maxSize;
  }
  else {        // multi-freedom case
    int theSkip = nSkips[theFreedom];
    int theBound = (sizes[theFreedom] - 1)*theSkip;
    int theDim = partDims[theFreedom];
    int k;
    // sum over other indices
    for(k=sizes[theFreedom],l=theBound; sum < epsilon; k--,l-=theSkip)
      for(i=0; i<theSkip; i++)
        for(j=0; j<maxSize; j+=theDim) {
          sum += norm(data[i+j+l]);    // add up probability
        }
    k += (padSize + 1);
    if( k > nDims[theFreedom] ) {  // enough states allocated?
      std::cerr << "Warning!  Significant probabilty in top states!" << std::endl;
      std::cerr << "k = " << k << std::endl;
      k = nDims[theFreedom];
    }
    if( k > sizes[theFreedom] ) stretchFreedom(theFreedom,k);
    if( k < sizes[theFreedom] ) shrinkFreedom(theFreedom,k);
  }
}

void State::stretchFreedom(int theFreedom, int theSize)
  //
  // If the space being used for a degree of freedom is not sufficient to
  // avoid truncation errors, stretch the available memory.
  //
{
  int i,j,k;
  int newDim,newSize;

  // Stretch the data into a larger space
  newDim = nSkips[theFreedom]*theSize;    // new partial dimension
  newSize = theSize*maxSize/sizes[theFreedom];  // new total size
  for(i=maxSize-partDims[theFreedom],j=newSize-newDim; i>=0;
      i-=partDims[theFreedom],j-=newDim) {
    for(k=partDims[theFreedom]-1; k>=0; k--)
      data[j+k] = data[i+k];      // stretch data
    for(k=newDim-1; k>=partDims[theFreedom]; k--)
      data[j+k] = 0;        // fill empty spaces w/zero
  }
  // Recompute nSkips and partDims
  sizes[theFreedom] = theSize;      // size of freedom
  partDims[theFreedom] = newDim;
  maxSize = newSize;        // total size used
  for(i=theFreedom+1; i<nFreedoms; i++) {
    nSkips[i] = partDims[i-1];      // new nSkips
    partDims[i] = nSkips[i]*sizes[i];    // new partDims
  }
  mySize = maxSize;
}

void State::shrinkFreedom(int theFreedom, int theSize)
  //
  // If more memory is being used than is necessary to avoid truncation
  // error, shrink memory used by freedom.
  //
{
  int i,j,k;
  int newDim,newSize;

  // Compact the data into a smaller space
  newDim = nSkips[theFreedom]*theSize;    // new partial dimension
  newSize = theSize*maxSize/sizes[theFreedom];  // new total size
  for(i=0,j=0; i<maxSize; i+=partDims[theFreedom],j+=newDim) //old and new outer indices
    for(k=0; k<newDim; k++)      // inner indices
      data[j+k] = data[i+k];      // compact data
  // Pad out the empty space with zeros
  for(i=newSize; i<maxSize; i++) data[i] = 0;
  // Recompute nSkips and partDims
  sizes[theFreedom] = theSize;      // size of freedom
  partDims[theFreedom] = newDim;
  maxSize = newSize;        // total size used
  for(i=theFreedom+1; i<nFreedoms; i++) {
    nSkips[i] = partDims[i-1];      // new nSkips
    partDims[i] = nSkips[i]*sizes[i];    // new partDims
  }
  mySize = maxSize;
}

void State::diagnostic()  // for debugging
{
  int i;
  std::cout << "Number of Freedoms " << nFreedoms << ".\n";
  std::cout << "Total Dimensions " << totalDim << ".\n";
  std::cout << "MaxSize " << maxSize << ".\n";
  std::cout << "My Type " << myType << ".\n";
  std::cout << "My Size " << mySize << ".\n";
  std::cout << "nSkip " << nSkip << ".\n";
  std::cout << "coord " << coord << ".\n";
  std::cout << "data pointer " << data << ".\n";
  std::cout << "My Pointer " << myPointer << ".\n";
  std::cout << "Dimensions:\n";
  if (nFreedoms > 1 ) {
    for (i=0; i<nFreedoms; i++)
      std::cout << "  " << nDims[i] << ".\n";
    std::cout << "Sizes:\n";
    for (i=0; i<nFreedoms; i++)
      std::cout << "  " << sizes[i] << ".\n";
    std::cout << "Skips:\n";
    for (i=0; i<nFreedoms; i++)
      std::cout << "  " << nSkips[i] << ".\n";
    std::cout << "Partial dimensions:\n";
    for (i=0; i<nFreedoms; i++)
      std::cout << "  " << partDims[i] << ".\n";
    std::cout << "Centers of coordinates:\n";
    for (i=0; i<nFreedoms; i++)
      std::cout << "  " << coords[i] << ".\n";
  }
}

// Friend I/O functions

std::ostream& operator<<( std::ostream& s, const State& a )
  //
  // Output a State to an output stream
  //
{
  int i;

  s << a.nFreedoms << std::endl;
  s << a.totalDim << std::endl;
  s << a.maxSize << std::endl;
  s << a.mySize << std::endl;
  s << a.nSkip << std::endl;
  s << a.coord << std::endl;
  s << a.myType << std::endl;
  if( a.nFreedoms > 1 )
    for(i=0; i<a.nFreedoms; i++) {
      s << a.nDims[i] << std::endl;
      s << a.sizes[i] << std::endl;
      s << a.nSkips[i] << std::endl;
      s << a.partDims[i] << std::endl;
      s << a.freedomTypes[i] << std::endl;
      s << a.coords[i] << std::endl;
    }
  for (i=0; i<a.maxSize; i++)
    s << a.data[i] << std::endl;
  return s;
}

std::istream& operator>>( std::istream& s, State& a )
  //
  // Input a State from an input stream
  //
{
  int i;

  a.free();      // free up State
  s >> a.nFreedoms;
  s >> a.totalDim;
  s >> a.maxSize;
  s >> a.mySize;
  s >> a.nSkip;
  s >> a.coord;
  s >> a.myType;
  if( a.totalDim > 0 ) {
    a.data = new Complex[a.totalDim];
    a.myPointer = a.data;
  }
  if( a.nFreedoms > 1 ) {
    a.nDims = new int[a.nFreedoms];
    a.sizes = new int[a.nFreedoms];
    a.nSkips = new int[a.nFreedoms];
    a.partDims = new int[a.nFreedoms];
    a.freedomTypes = new FreedomType[a.nFreedoms];
    a.coords = new Complex[a.nFreedoms];
    for(i=0; i<a.nFreedoms; i++) {
      s >> a.nDims[i];
      s >> a.sizes[i];
      s >> a.nSkips[i];
      s >> a.partDims[i];
      s >> a.freedomTypes[i];
      s >> a.coords[i];
    }
  }
  for (i=0; i<a.maxSize; i++)
    s >> a.data[i];
  return s;
}

// apply a PrimaryOperator to the State

void State::apply(PrimaryOperator& theOp, int hc, int theFreedom, FreedomType theType, double t)
{
#ifdef DEBUG_TRACE
  std::cout << "State::apply(theOp, hc="<<hc<<", theFreedom="<<theFreedom<<", FreedomType="<<theType<<", t="<<t<<") entered" << std::endl;
#endif
  int i,j;
  if( nFreedoms == 1 ) {    // 1 freedom case
#ifndef OPTIMIZE_QSD
    if( (theType != myType) && (theType != ALL) )
      error("Incompatible operator type in State::apply!");
#endif
    theOp.applyTo(*this,hc,t);
  }
  else {        // multi-freedom case
#ifndef OPTIMIZE_QSD
    if( (theFreedom < 0) || (theFreedom >= nFreedoms) )
      error("Illegal degree of freedom requested in State::apply!");
    if( (theType != freedomTypes[theFreedom]) && (theType != ALL) )
      error("Incompatible operator type in State::apply!");
#endif
    mySize = sizes[theFreedom];    // set up 1D state for PrimaryOperator
    nSkip = nSkips[theFreedom];
    coord = coords[theFreedom];
    for (i=0; i<nSkip; i++) {        // loop over inner and
      for (j=0; j<maxSize; j += partDims[theFreedom]) {  // outer indices
        myPointer = data + i + j;
        if( i+j > maxSize ) error("Gone off end of array!");
        theOp.applyTo(*this,hc,t);
      } }
    mySize = maxSize;
    nSkip = 1;
    coord = coords[0];
    myPointer = data;
  }
}

State operator*(const Complex& a, const State& b) {  // multiplication
  State result(b);
  return (result *= a);
}

State operator*(const State& b, const Complex& a) {  // multiplication
  State result(b);
  return (result *= a);
}

State operator*(double a, const State& b) {    // multiplication
  State result(b);
  return (result *= a);
}

State operator*(const State& b, double a) {    // multiplication
  State result(b);
  return (result *= a);
}

State operator*(ImaginaryUnit a, const State& b) {  // multiplication
  State result(b);
  return (result *= a);
}

State operator*(const State& b, ImaginaryUnit a) {  // multiplication
  State result(b);
  return (result *= a);
}

State operator+(const State& a, const State& b) {  // addition
  State result(a);
  return (result += b);
}

State operator-(const State& a, const State& b) {  // subtraction
  State result(a);
  return (result -= b);
}

State operator+(const State& a) {      // unary +
  State result(a);
  return result;
}

State operator-(const State& a) {      // unary -
  State result(a);
  for (int i=0; i<result.maxSize; i++) result.data[i] = - result.data[i];
  return result;
}

int State::size() { return mySize; }  // return state size

int State::getSize(int theFreedom)
{
#ifndef OPTIMIZE_QSD
  if( (theFreedom < 0) || (theFreedom >= nFreedoms) )
    error("Illegal freedom passed to State::getSize!");
#endif
  return sizes[theFreedom];
}

Complex State::centerCoords() { return coord; }  // return center of coordinates

Complex State::getCoords(int theFreedom)
{
#ifndef OPTIMIZE_QSD
  if( (theFreedom<0) || (theFreedom >= nFreedoms) )
    error("Illegal freedom requested in State::getCoords!");
#endif
  if( nFreedoms > 1 ) return coords[theFreedom];
  return coord;
}

void State::setCoords(Complex& alpha, int theFreedom)
  //
  // Reset the value of the phase space coordinate for a freedom.
{
#ifndef OPTIMIZE_QSD
  if( (theFreedom<0) || (theFreedom >= nFreedoms) )
    error("Illegal freedom requested in State::setCoords!");
#endif
  if( nFreedoms == 1 )
    coord = alpha;
  else {
    coords[theFreedom] = alpha;
    if( theFreedom == 0 ) coord = alpha;
  }
}

void State::displaceCoords(Complex& alpha, int theFreedom)
  //
  // Reset the value of the phase space coordinate for a freedom.
{
#ifndef OPTIMIZE_QSD
  if( (theFreedom<0) || (theFreedom >= nFreedoms) )
    error("Illegal freedom requested in State::setCoords!");
#endif
  if( nFreedoms == 1 )
    coord += alpha;
  else {
    coords[theFreedom] += alpha;
    if( theFreedom == 0 ) coord += alpha;
  }
}

// Move the center of coordinates to a new position in phase space.

void State::moveToCoords(const Complex& alpha, int theFreedom, double shiftAccuracy)
{
#ifndef OPTIMIZE_QSD
  if( (theFreedom < 0) || (theFreedom >= nFreedoms) )
    error("Illegal degree of freedom passed in State::moveToCoords!");
  if( nFreedoms == 1 ) {
    if( myType != FIELD )
      error("Can only change basis of FIELD freedoms in State::moveToCoords!");
  }
  else {
    if( freedomTypes[theFreedom] != FIELD )
      error("Can only change basis of FIELD freedoms in State::moveToCoords!");
  }
#endif
  Complex displacement;
  displacement = getCoords(theFreedom) - alpha;
  moveCoords(displacement,theFreedom,shiftAccuracy);
}

// Move the center of coordinates to match those of another state.

void State::centerOn(State& psi, double shiftAccuracy)
{
#ifndef OPTIMIZE_QSD
  if( nFreedoms != psi.nFreedoms )
    error("Mismatched number of degrees of freedom in State::centerOn!");
#endif
  Complex alpha;
  if( nFreedoms == 1 ) {
    if( (myType != FIELD) || (psi.myType != FIELD) )
      error("Can't change basis of non-FIELD type in State::centerOn!");
    alpha = centerCoords() - psi.centerCoords();
    moveCoords(alpha,0,shiftAccuracy);
  }
  else {
    for (int i=0; i<nFreedoms; i++) {
      if( freedomTypes[i] != psi.freedomTypes[i] )
        error("Mismatched types in State::centerOn!");
      if( freedomTypes[i] == FIELD ) {
        alpha = getCoords(i) - psi.getCoords(i);
        moveCoords(alpha,i);
      }
    }
  }
}

// Recenter the coordinates at the current mean position in phase space.

void State::recenter(int theFreedom, double shiftAccuracy)
{
#ifndef OPTIMIZE_QSD
  if( (theFreedom >= nFreedoms) || (theFreedom < 0) )
    error("Illegal degree of freedom in State::recenter!");
  if( nFreedoms == 1 ) {
    if (myType != FIELD) error("Only type FIELD can be recentered!");
  }
  else {
    if (freedomTypes[theFreedom] != FIELD)
      error("Only type FIELD can be recentered!");
  }
#endif
  int i,j,k,l;
  State dummy=(*this);
  Complex magnitude, alpha;
  magnitude = (*this)*(*this);
  if( magnitude == 0 ) error("Can't recenter zero state!");
  if( nFreedoms == 1 ) {    // 1 freedom case
    if( mySize > intSqrtMax ) setIntSqrtMax(mySize);
    for (k=1; k<mySize; k++)  // apply local annihilation op.
      dummy.myPointer[k-1] = intSqrtList[k]*dummy.myPointer[k];
    dummy.myPointer[mySize-1] = 0;
  }
  else {        // multi-freedom case
    mySize = sizes[theFreedom];
    nSkip = nSkips[theFreedom];
    if( mySize > intSqrtMax ) setIntSqrtMax(mySize);
    for (i=0; i<nSkip; i++)    // apply local annihilation op
      for (j=0; j<maxSize; j+=partDims[theFreedom]) {  // loop over other
        dummy.myPointer = dummy.data + i + j;    // indices
        for (k=1,l=nSkip; k<mySize; k++,l+=nSkip)
          dummy.myPointer[l-nSkip] = intSqrtList[k]*dummy.myPointer[l];
        dummy.myPointer[(mySize-1)*nSkip] = 0;
      }
    mySize = maxSize;
    nSkip = 1;
    dummy.myPointer = dummy.data;
  }
  alpha = - ((*this)*dummy)/magnitude;  // get relative coords in phase space
  if( abs(alpha) > shiftAccuracy )
    moveCoords(alpha,theFreedom, shiftAccuracy);  // shift basis
}

// Move the center of coordinates by a displacement alpha in
// phase space.

void State::moveCoords(const Complex& alpha, int theFreedom, double shiftAccuracy)
{
#ifndef OPTIMIZE_QSD
  if( (theFreedom >= nFreedoms) || (theFreedom < 0) )
    error("Illegal freedom passed to State::moveCoords!");
#endif
  Complex mag=(*this)*(*this);
  if( mag == 0 ) {        // move zero state
    if( nFreedoms == 1 ) coord -= alpha;
    if( nFreedoms > 1 ) {
      coords[theFreedom] -= alpha;
      coord = coords[0];
    }
    return;
  }
  int numshifts;        // Break basis shift into
  numshifts = int(abs(alpha)/shiftAccuracy)+1;  // sequence of small shifts
  Complex delta=alpha/numshifts;    // according to accuracy
  Complex deltacc=conj(delta);      // requirement.
  Complex theCoord, phase=0;
  if (nFreedoms == 1) {        // 1 freedom case
#ifndef OPTIMIZE_QSD
    if( myType != FIELD) error("Can only move coordinates for FIELD type.");
#endif
    theCoord = coord;
    if( mySize > intSqrtMax) setIntSqrtMax(mySize);
    for (int i=0; i<numshifts; i++) {
      moveStep(delta,deltacc);      // 1 small shift
      phase += imag(delta*conj(theCoord));  // sum up phase changes
      theCoord -= delta;
    }
    coord = theCoord;
  }
  else {          // multi-freedom case
#ifndef OPTIMIZE_QSD
    if( freedomTypes[theFreedom] != FIELD )
      error("Can only move coordinates for FIELD type.");
#endif
    mySize = sizes[theFreedom];
    nSkip = nSkips[theFreedom];
    theCoord = coords[theFreedom];
    if( mySize > intSqrtMax) setIntSqrtMax(mySize);
    for (int i=0; i<numshifts; i++) {
      for (int j=0; j<nSkip; j++)    // loop over other indices
        for (int k=0; k<maxSize; k+=partDims[theFreedom]) {
          myPointer = data + j + k;
          moveStep(delta,deltacc);    // 1 small shift
        }
      phase += imag(delta*conj(theCoord));  // sum up phase changes
      theCoord -= delta;
    }
    coords[theFreedom] = theCoord;
    if( theFreedom == 0 ) coord = coords[0];
    myPointer = data;
    mySize = maxSize;
    nSkip = 1;
  }
  phase.timesI();
  (*this) *= exp(phase);      // remove phase change
  normalize();          // normalize state
}

// Shift the center of coordinates by a single infinitesimal step

void State::moveStep(Complex& delta, Complex& deltacc)
{
  int i,j;
  Complex psizero,psinth,temp,previous;
  psizero = (*this)[0] - deltacc*(*this)[1];
  psinth = (*this)[mySize-1] + delta*intSqrtList[mySize-1]*(*this)[mySize-2];
  previous = (*this)[0];
  for (i=1,j=nSkip; i<(mySize-1); i++,j+=nSkip) {  // infinitesimal
    temp = myPointer[j];        // displacement
    myPointer[j] += delta*intSqrtList[i]*previous;
    myPointer[j] -= deltacc*intSqrtList[i+1]*myPointer[j+nSkip];
    previous = temp;
  }
  (*this)[0] = psizero;
  (*this)[mySize-1] = psinth;
}

void State::error(const char* message) const  // error handling
{
  std::cerr << message << "\n";
  exit(1);
}

