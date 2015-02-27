// onespin.cc

// One spin, one or several trajectories.

#include <math.h>
#include <stdio.h>
#include <iostream.h>
#include <fstream.h>

#include "ACG.h"
#include "Traject.h"
#include "State.h"
#include "Operator.h"
#include "SpinOp.h"
#include "Complex.h"

main()
{
    // Basic operators
  IdentityOperator id;
  NullOperator null;
  SigmaX sx;
  SigmaY sy;
  SigmaZ sz;
  SigmaPlus sp;
  Operator sm = sp.hc();    // sp.hc() is the Hermitian conjugate of sp.

    // The Hamiltonian
  double omega=0.5;
  double epsilon=0.1;
  Operator H = omega*sz + epsilon*sx;

    // The Lindblad operators
  const int nOfLindblads = 1;
  double gamma=0.1;
  Operator L1 = gamma*sm;
  Operator L[nOfLindblads] = {L1};

    // The initial state
  State psi0(2,SPIN);			// ground state ("spin down")
  psi0 *= sp;                           // excited state ("spin up")

    // The random number generator
  int seed = 74298;                 // change seed for independent runs
  ACG gen(seed,55);                 // don't change the value 55
  ComplexNormal rand1(&gen);

    // Stepsize and integration time  
  double dt=0.01;   // basic time step
  int numdts=10;    // time interval between outputs = numdts*dt
  int numsteps=20;  // total integration time = numsteps*numdts*dt
  
  double accuracy = 0.000001;

  AdaptiveStep theStepper(psi0,H,nOfLindblads,L,accuracy);
    // deterministic part: adaptive stepsize 4th/5th order Runge Kutta
    // stochastic part: fixed stepsize Euler

//AdaptiveStochStep theStepper(psi0,H,nOfLindblads,L,accuracy);
    // deterministic part: adaptive stepsize 4th/5th order Runge Kutta
    // stochastic part: Euler, same (variable) stepsize as deterministic part

    // Output
  const int nOfOut = 2;
  Operator outlist[nOfOut] = {sz,sx};	        // Operators to output

  char *flist[nOfOut] = {"sz.out","sx.out"};    // Output files
    // While the program is running, in addition to the data written in the
    // output files, 7 columns are written to standard output:
    //  the time `t' in column 1;
    //  4 values determined by the array `pipe' (see below) in columns 2-5;
    //  the effective dimension of Hilbert space in column 6;
    //  the number of adaptive steps taken in column 7. 

  int pipe[4] = {1,3,5,7};
    // The 4 numbers in `pipe' refer to a list formed by the real 
    // and imaginary parts of the expectations and variances of the operators
    // in `outlist'. 
    // Example: For outlist={sz,sx} as above, pipe={1,3,5,7} refers to 
    // the 1st, 3rd, 5th and 7th entries in the list
    //  { Re(<sz>),Im(<sz>),Re(<sz^2>-<sz>^2),Im(<sz^2>-<sz>^2),
    //     Re(<sx>),Im(<sx>),Re(<sx^2>-<sx>^2),Im(<sx^2>-<sx>^2) }, 
    // i.e. to the values
    //  Re(<sz>), Re(<sz^2>-<sz>^2), Re(<sx>), Re (<sx^2>-<sx>^2).

    // Integrate `nTraj' trajectories, all starting from the same 
    //   initial state `psi0'.
  int nTraj = 1;
  for( int i=0; i<nTraj; i++ ) {
    Trajectory theTraject(psi0,dt,theStepper,&rand1);
    theTraject.plotExp(nOfOut,outlist,flist,pipe,numdts,numsteps);
    cout << endl;
  }
}






