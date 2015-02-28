// Simple harmonic oscillator, using the moving basis, one trajectory.
// Different unravelings can be chosen by uncommenting or
// commenting out the appropriate lines.

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>

#include "ACG.h"
#include "Traject.h"
#include "State.h"
#include "Operator.h"
#include "FieldOp.h"
#include "Complex.h"

main()
{
    // Basic operators
  IdentityOperator Id;
  AnnihilationOperator A;
  NumberOperator N;
  Operator Ac = A.hc();

    // The Hamiltonian
  double omega = 1.0;
  Operator H = omega*Ac*A;

    // The Lindblad operators
  double gamma = 1.0;
  const int nOfLindblads = 1;
  Operator L1 = gamma*A;
  Operator L[nOfLindblads] = {L1};

    // The initial state
  int dim = 4;           // cutoff Hilbert space dimension, can be small
                         //  since the moving basis is used.
  Complex alpha(1.4,-0.4);

//int m=5;
//State psi(dim,m,alpha); // excited coherent state |alpha,m>

  State psi(dim,0,alpha); // coherent state |alpha> 
                          // (centered at alpha, i.e., represented as the
                          // ground state of the excited basis {|alpha,m>})
                          // See "State.h" for more details.

    // The random number generator
  int seed = 74298;
  ACG gen(seed,55);
  ComplexNormal rand1(&gen);   // noise function for QSD
//ComplexUniform rand1(&gen);  // noise function for Quantum Jumps

    // Stepsize and integration time  
  double dt=0.01;   // basic time step
  int numdts=10;    // time interval between outputs = numdts*dt
  int numsteps=100; // total integration time = numsteps*numdts*dt
  
  double accuracy = 0.000001;

  AdaptiveStep theStepper(psi,H,nOfLindblads,L,accuracy);
    // QSD unraveling
    // deterministic part: adaptive stepsize 4th/5th order Runge Kutta
    // stochastic part: fixed stepsize Euler

//AdaptiveStochStep theStepper(psi,H,nOfLindblads,L,accuracy);
    // QSD unraveling
    // deterministic part: adaptive stepsize 4th/5th order Runge Kutta
    // stochastic part: Euler, 
    //     same (variable) stepsize as deterministic part

//AdaptiveJump theStepper(psi,H,nOfLindblads,L,accuracy,"jump_times");
    // Quantum Jumps unraveling (jumps in file "jump_times")
    // deterministic part: adaptive stepsize 4th/5th order Runge Kutta
    // stochastic part: Euler, 
    //      same (variable) stepsize as deterministic part

  int nOfMovingFreedoms = 1;   // Moving basis for 1 degree of freedom.

// We dynamically adjust the number of basis vectors. Our criterion for
// this adjustment depends on parameters `pCutoff' the cutoff
// probability, and `nPad', the pad size, which represents the number of
// boundary basis states that are checked for significant probability. We
// require the total probability of the top `nPad' states to be no
// greater than `pCutoff', increasing and decreasing the number of states
// actually used accordingly, as the integration proceeds along the
// quantum trajectory.
// See  J. Phys. A 28, 5401 (1995) for more details.

  double pCutoff = 0.01;
  int nPad = 2;

    // Output
  const int nOfOut = 2;
  Operator outlist[nOfOut] = {A,N};	    // Operators to output

  char *flist[nOfOut] = {"A.out","N.out"};  // Output files

  int pipe[4] = {1,2,5,7};
    // Standard output:
    // t,Re<A>,Im<A>,<N>,<N^2>-<N>^2,dim,steps
    //  where `t' is time, 
    //  `dim' is the effective dimension of Hilbert space,
    //  and steps is the number of adaptive steps taken.
    // (for more explanation see `onespin.cc')

    // Simulate one trajectory (for several trajectories, see `onespin.cc')
  Trajectory theTraject(psi,dt,theStepper,&rand1);
  theTraject.plotExp(nOfOut,outlist,flist,pipe,numdts,numsteps,
		     nOfMovingFreedoms, pCutoff, nPad);
}
