// Uses many features, may serve as a template.

#include "Complex.h"
#include "ACG.h"
#include "CmplxRan.h"
#include "State.h"
#include "Operator.h"
#include "FieldOp.h"
#include "SpinOp.h"
#include "Traject.h"

int main() 
{
// Primary Operators
  AnnihilationOperator A1(0);  // 1st freedom
  NumberOperator N1(0);
  IdentityOperator Id1(0);
  AnnihilationOperator A2(1);  // 2nd freedom
  NumberOperator N2(1);
  IdentityOperator Id2(1);
  SigmaPlus Sp(2);             // 3rd freedom
  IdentityOperator Id3(2);
  Operator Sm = Sp.hc();       // Hermitian conjugate
  Operator Ac1 = A1.hc();
  Operator Ac2 = A2.hc();
// Hamiltonian
  double E = 20.0;           
  double chi = 0.4;      
  double omega = -0.7;       
  double eta = 0.001;
  Complex I(0.0,1.0);
  Operator H = (E*I)*(Ac1-A1)
             + (0.5*chi*I)*(Ac1*Ac1*A2 - A1*A1*Ac2)
             + omega*Sp*Sm + (eta*I)*(A2*Sp-Ac2*Sm);
// Lindblad operators
  double gamma1 = 1.0;       
  double gamma2 = 1.0;       
  double kappa = 0.1;        
  const int nL = 2;
  Operator L[nL]={sqrt(2*gamma1)*A1,sqrt(2*gamma2)*A2};
// Initial state
  State phi1(50,FIELD);       // see paper Section 4.2
  State phi2(50,FIELD);
  State phi3(2,SPIN);
  State stateList[2] = {phi1,phi2};
  State psiIni(2,stateList);
// Trajectory
  double dt = 0.01;    // basic time step                            
  int numdts = 100;    // time interval between outputs = numdts*dt  
  int numsteps = 5;    // total integration time = numsteps*numdts*dt
  int nOfMovingFreedoms = 2;
  double epsilon = 0.01;     // cutoff probability
  int nPad = 2;              // pad size
  ACG gen(38388389);         // random number generator with seed
  ComplexNormal rndm(&gen);  // Complex Gaussian random numbers
  AdaptiveStep stepper(psiIni, H, nL, L);       // see paper Section 5
// Output
  const int nOfOut = 2;
  Operator outlist[nOfOut]={ Ac1*A1, Ac2*A2 };
  char *flist[nOfOut]={"X1.out","X2.out"};
  int pipe[] = {1,2,5,7}; // controls standard output (see `onespin.cc')
// Simulate one trajectory (for several trajectories see `onespin.cc')
  Trajectory traj(psiIni, dt, stepper, &rndm);  // see paper Section 5
  traj.plotExp( nOfOut, outlist, flist, pipe, numdts, numsteps,
                nOfMovingFreedoms, epsilon, nPad );
}

