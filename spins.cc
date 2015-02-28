// Two interacting spins, time-dependent Hamiltonian, one trajectory.

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>

#include "ACG.h"
#include "Traject.h"
#include "State.h"
#include "Operator.h"
#include "SpinOp.h"
#include "Complex.h"

#define TWOPI 6.2831853

double driving(double t)           // Step function.
{
  if ((t<1) && (t>0)) return 1.0;
  return 0.0;
}

main()
{
    // Basic operators
  IdentityOperator id1(0);  // Operates on the zeroth degree of freedom.
  IdentityOperator id2(1);  // Operates on the first degree of freedom.
  SigmaX sx1(0);
  SigmaX sx2(1);
  SigmaY sy1(0);
  SigmaY sy2(1);
  SigmaZ sz1(0);
  SigmaZ sz2(1);
  SigmaPlus sp1(0);
  SigmaPlus sp2(1);
  Operator sm1 = sp1.hc();  // sp1.hc() is the Hermitian conjugate of sp1.
  Operator sm2 = sp2.hc();

    // The Hamiltonian
  RealFunction f = driving;
  double OMEGA = 0;
  double KAPPA=TWOPI/2;
  Operator H = OMEGA*sz1 + OMEGA*sz2 + KAPPA*(f*(sm1*sp2 + sm2*sp1));

    // The Lindblad operators
  double GAMMA=0.3;
  const int nOfLindblads = 2;
  Operator L1 = GAMMA*sm1;
  Operator L2 = GAMMA*sm2;
  Operator L[nOfLindblads] = {L1,L2};

    // The initial state
  State psi1(2,SPIN);			// 1-freedom state ("down")
  State psi2(2,SPIN);			// 1-freedom state ("down")
  State psilist[2] = {psi1,psi2};       
  State psi(2,psilist);                 // Product state  ("down,down")
  psi *= sp2;         // similar to `psi=sp2*psi', results in ("down,up") 
    // To create an entangled state, add two or more product states.

    // The random number generator
  int seed = 74298;
  ACG gen(seed,55);
  ComplexNormal rand1(&gen);

    // Stepsize and integration time  
  double dt=0.01;   // basic time step
  int numdts=10;    // time interval between outputs = numdts*dt
  int numsteps=20;  // total integration time = numsteps*numdts*dt
  
  double accuracy = 0.000001;

  AdaptiveStep theStepper(psi,H,nOfLindblads,L,accuracy);
    // deterministic part: adaptive stepsize 4th/5th order Runge Kutta
    // stochastic part: fixed stepsize Euler

    // Output
  const int nOfOut = 2;
  Operator outlist[nOfOut] = {sz1,sz2};	        // Operators to output

  char *flist[nOfOut] = {"sz1.out","sz2.out"};  // Output files

  int pipe[4] = {1,3,5,7};
    // Controls standard output. For details see `onespin.cc`.

    // Integrate one trajectory (several trajectories see `onespin.cc')
  Trajectory theTraject(psi,dt,theStepper,&rand1);
  theTraject.plotExp(nOfOut,outlist,flist,pipe,numdts,numsteps);
}
