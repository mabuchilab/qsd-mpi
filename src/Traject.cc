//   Traject.cc	-- Stochastic simulation of QSD trajectories.
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream.h>
#include <fstream.h>

#include "Complex.h"
#include "State.h"
#include "Operator.h"
#include "ACG.h"
#include "Normal.h"
#include "Traject.h"

#define MAXSTP 10000
#define TINY 1.0e-30
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

static const char rcsid[] = "$Id: Traject.cc,v 3.5 1996/11/22 10:34:37 mrigo Exp $";

Trajectory::Trajectory()	// Constructor
{
  error("A trajectory must be constructed with initial data!");
}

AdaptiveStep::AdaptiveStep(const State& psi, const Operator& theH, int theNL,
	const Operator* theL, double theEpsilon)
//
// Constructor for AdaptiveStep; this contains all the information needed
// to integrate the QSD equation for a single adaptive timestep.  The
// Hamiltonian and Lindblad operators are stored as private data; the
// State psi is used to initialize the temporary States used in the
// Runge-Kutta integration; nL is the number of Lindblad operators, and
// is also stored as private data; and the integration is performed up to
// an accuracy epsilon.
//
{
  H = theH;				// store private data
  nL = theNL;
  L = new Operator[nL];
  Ldag = new Operator[nL];
  dxi = new Complex[nL];
  for (int i=0; i<nL; i++) {
    L[i] = theL[i];
    Ldag[i] = L[i].hc();
    dxi[i] = 0;
  }
  epsilon = theEpsilon;
  stochasticFlag=0;
  numdtsUsed = 0;
  dtNumList = 0;
  listDim = 0;
  listPtr = 0;

  ak2 = ak3 = ak4 = ak5 = psi;		// initialize temporaries
  ytemp = yerr = y = yout = dydt = psi;
  temp1 = newsum = psi;
}

AdaptiveStep::~AdaptiveStep()	// AdaptiveStep destructor
{
  if( listDim != 0 )
#ifndef NON_GNU_DELETE
    delete[] dtNumList;		// deallocate memory
#else
    delete[listDim] dtNumList;	// deallocate memory
#endif
}

AdaptiveJump::~AdaptiveJump()	// AdaptiveJump destructor
{
  if( nL != 0 )
#ifndef NON_GNU_DELETE
    delete[] lindVal;		// deallocate memory
#else
    delete[nL] lindVal;		// deallocate memory
#endif
  if( outFile != 0 ) fclose(outFile);
}

AdaptiveOrthoJump::~AdaptiveOrthoJump()     // AdaptiveOrthoJump destructor
{
  if( nL != 0 )
#ifndef NON_GNU_DELETE
    delete[] lindVal;		// deallocate memory
#else
    delete[nL] lindVal;		// deallocate memory
#endif
  if( outFile != 0 ) fclose(outFile);
}

Trajectory::Trajectory(const State& thePsiIni, double thedt,
  IntegrationStep& theStepper, ComplexRandom* theRand, double theT0)
//
// Constructor for trajectory class.  This stores all the
// information needed to calculate a QSD trajectory.  The initial state,
// size of the output step, and number of dts per output step are stored
// as member data, and the size of dt is calculated; the seed for the
// random number generator is stored, as is the initial time; and the
// IntegrationStep contains all the information about the Hamiltonian and
// Lindblad operators, as well as the integration algorithm.
//
{
  psi = thePsiIni;
  dt = thedt;
  rndm = theRand;
  t0 = theT0;
  stepper = &theStepper;
}

IntegrationStep::~IntegrationStep()	// destructor
{
#ifndef NON_GNU_DELETE
  delete [] L;				// dispose of allocated memory
  delete [] Ldag;
  delete [] dxi;
#else // ---- NON GNU CODE ----
  delete [nL] L;			// dispose of allocated memory
  delete [nL] Ldag;
  delete [nL] dxi;
#endif // ---- NON GNU CODE ----
}

void AdaptiveStep::operator()(State& psi, double t, double dt,
	double& dtlast, ComplexRandom* rndm)
//
// This member function calculates a single QSD integration step.  It
// invokes the member function odeint (adapted from the Numerical Recipes
// routine) to advance the deterministic part of the equation in time,
// and uses the supplied function rndm to advance the stochastic
// part (with noise dxi).  The state psi is evolved from t to t+dt; dtlast
// contains the actual physical stepsize used by the adaptive stepsize
// routine.  (In case of a constant stepsize stepper, dtlast=dt.)
//
{
  int nok, nbad;

  if( rndm != 0 ) {				// any noise?
    double sqrtdt = sqrt(dt);
    stochasticFlag = 1;
	// set flag for derivs to generate stoch. part and store in newsum
    for (int i=0; i<nL; i++) {			// for each Lindblad...
      dxi[i] = sqrtdt*(*rndm)();		// ...independent noise
    }
  }
  newsum = psi;					// initialize temp
  newsum = 0;
  odeint(psi,t,t+dt,epsilon,dtlast,0.0,nok,nbad,dxi);
	// generate deterministic evolution (and also stochastic term)
  psi += newsum;				// add stochastic term
  psi.normalize();
}

void AdaptiveStochStep::operator()(State& psi, double t, double dt,
	double& dtlast, ComplexRandom* rndm)
//
// This member function calculates a single QSD integration step.  It
// invokes the member function odeint (adapted from the Numerical Recipes
// routine) to advance the deterministic part of the equation in time,
// and uses the supplied function rndm to advance the stochastic
// part (with noise dxi).  The state psi is evolved from t to t+dt; dtlast
// contains the actual physical stepsize used by the adaptive stepsize
// routine.  (In case of a constant stepsize stepper, dtlast=dt.)
//
{
  int nok, nbad;
  cout << " Using adaptiveStochStep Rad " << endl;
  newsum = psi;					// initialize temp
  odeint(psi,t,t+dt,epsilon,dtlast,0.0,nok,nbad,rndm);
	// generate deterministic evolution (and also stochastic)
}

void AdaptiveJump::operator()(State& psi, double t, double dt,
	double& dtlast, ComplexRandom* rndm)
//
// This member function calculates a single Quantum Jumps integration step.
// It invokes the member function odeint (adapted from the Numerical Recipes
// routine) to advance the deterministic part of the equation in time,
// and uses the supplied function rndm to advance the stochastic
// part (with noise dxi).  The state psi is evolved from t to t+dt; dtlast
// contains the actual physical stepsize used by the adaptive stepsize
// routine.  (In case of a constant stepsize stepper, dtlast=dt.)
//
{
  int nok, nbad;

  odeint(psi,t,t+dt,epsilon,dtlast,0.0,nok,nbad,rndm);
	// generate deterministic evolution (and also stochastic)
}

void AdaptiveOrthoJump::operator()(State& psi, double t, double dt,
	double& dtlast, ComplexRandom* rndm)
//
// This member function calculates a single QSD Quantum Jumps integration step.
// QSD Quantum Jumps as the same deterministic part as QSD but the stochastic
// part correspond to jumps.
// It invokes the member function odeint (adapted from the Numerical Recipes
// routine) to advance the deterministic part of the equation in time,
// and uses the supplied function rndm to advance the stochastic
// part (with noise dxi).  The state psi is evolved from t to t+dt; dtlast
// contains the actual physical stepsize used by the adaptive stepsize
// routine.  (In case of a constant stepsize stepper, dtlast=dt.)
//
{
  int nok, nbad;

  odeint(psi,t,t+dt,epsilon,dtlast,0.0,nok,nbad,rndm);
	// generate deterministic evolution (and also stochastic)
}

State Trajectory::getState()
// Returns the value of psi stored in the trajectory.
{
  return psi;
}

void Trajectory::plotExp_obsolete( int nX, const Operator* X, FILE** fp, int* pipe,
	int dtsPerStep, int nOfSteps, int move,
	double delta, int width, double moveEps)
//
// This member function calculates a trajectory consisting of nOfSteps
// output steps, and at each step stores the expectation values of
// the nX operators X into files.  The array `pipe'
// instructs plotExp which expectation values it should print to the
// standard output.  The variable move contains the number of freedoms
// for which the moving basis and changing cutoff algorithms should be used;
// delta is the probability threshold for the cutoff adjustment, and
// width is the size of the "pad" for the cutoff.  After completing the
// trajectory, the initial state psi is reset to the final state.
//
{
  double t=t0;				// time
  double dtlast=dt;			// keep track of numerical timestep
  double xpipe[4];			// output to standard output
  if( nX > 0 ) {
    for( int k=0; k<4; k++ )
      if( pipe[k] < 1 || pipe[k] > 4*nX )
        cerr << "Warning: illegal argument 'pipe' in Trajectory::plotExp." 
          << endl;
  }

  State psi1 = psi;			// temporary state

  for( int i=0; i<nX; i++ ) {		// compute expectation values for
    psi1 = psi;				// output...
    psi1 *= X[i];
    Complex expec = psi*psi1;
    psi1 *= X[i];
    Complex expec2 = psi*psi1 - expec*expec;
    fprintf( fp[i], "%lG %lG %lG %lG %lG\n", 
	    t, real(expec), imag(expec), real(expec2), imag(expec2) ); // file
    for( int k=0; k<4; k++ ) {
      if( pipe[k] == 4*i+1 ) xpipe[k] = real(expec);          // standard i/o
      else if( pipe[k] == 4*i+2 ) xpipe[k] = imag(expec);
      else if( pipe[k] == 4*i+3 ) xpipe[k] = real(expec2);
      else if( pipe[k] == 4*i+4 ) xpipe[k] = imag(expec2);
    }
  }
  if( nX > 0 )
    printf("%lG %lG %lG %lG %lG %d %d\n",
      t,xpipe[0],xpipe[1],xpipe[2],xpipe[3],psi.size(),(*stepper).getNumdts());

  for (int n=1; n<=nOfSteps; n++) {
    for (int n1=0; n1<dtsPerStep; n1++) {
      (*stepper)(psi,t,dt,dtlast,rndm);
      t += dt;
      for(int k=0; k<move; k++) {		// use moving basis & cutoff
        psi.adjustCutoff(k,delta,width);
        psi.recenter(k,moveEps);
      }
    }
    for( int i=0; i<nX; i++ ) {		// compute expectation values for
      psi1 = psi;			// output...
      psi1 *= X[i];
      Complex expec = psi*psi1;
      psi1 *= X[i];
      Complex expec2 = psi*psi1 - expec*expec;
      fprintf( fp[i], "%lG %lG %lG %lG %lG\n", 
	      t, real(expec),imag(expec), real(expec2),imag(expec2) ); // file
      for( int k=0; k<4; k++ ) {
	if( pipe[k] == 4*i+1 ) xpipe[k] = real(expec);          // standard i/o
	else if( pipe[k] == 4*i+2 ) xpipe[k] = imag(expec);
	else if( pipe[k] == 4*i+3 ) xpipe[k] = real(expec2);
	else if( pipe[k] == 4*i+4 ) xpipe[k] = imag(expec2);
      }
    }
    if( nX > 0 )
      printf("%lG %lG %lG %lG %lG %d %d\n",
        t,xpipe[0],xpipe[1],xpipe[2],xpipe[3],psi.size(),(*stepper).getNumdts());
  }
//
// Save state to output file
//
//  ofstream os("saved_state.dat", ios::out );
//  if (!os) {
//    cerr << "Can't open output file:  saved_state.dat" << endl;
//    exit(1);
//  }
//  os << psi;
}


void Trajectory::plotExp( int nX, const Operator* X, char** fname, int* pipe,
	int dtsPerStep, int nOfSteps, int move,
	double delta, int width, double moveEps, char* savedState)
//
// This member function calculates a trajectory consisting of nOfSteps
// output steps, and at each step stores the expectation values of
// the nX operators X into files.  The array `pipe'
// instructs plotExp which expectation values it should print to the
// standard output.  The variable move contains the number of freedoms
// for which the moving basis and changing cutoff algorithms should be used;
// delta is the probability threshold for the cutoff adjustment, and
// width is the size of the "pad" for the cutoff.  After completing the
// trajectory, the initial state psi is reset to the final state.
// Unless the string variable `savedState' is equal to the Null pointer 
// (which is the default), the final state is saved to the file `savedState'.
{
  double t=t0;				// time
  double dtlast=dt;			// keep track of numerical timestep
  double xpipe[4];			// output to standard output
  FILE** fp;
  fp = new FILE*[nX];                   // nX files are used

  if( nX > 0 ) 
    for ( int i=0; i<nX; i++) fp[i] = fopen(fname[i],"w"); 

  if( nX > 0 ) {
    for( int k=0; k<4; k++ )
      if( pipe[k] < 1 || pipe[k] > 4*nX )
        cerr << "Warning: illegal argument 'pipe' in Trajectory::plotExp." 
          << endl;
  }

  State psi1 = psi;			// temporary state

  for( int i=0; i<nX; i++ ) {		// compute expectation values for
    psi1 = psi;				// output...
    psi1 *= X[i];
    Complex expec = psi*psi1;
    psi1 *= X[i];
    Complex expec2 = psi*psi1 - expec*expec;
    fprintf( fp[i], "%lG %lG %lG %lG %lG\n", 
	    t, real(expec), imag(expec), real(expec2), imag(expec2) ); // file
    for( int k=0; k<4; k++ ) {
      if( pipe[k] == 4*i+1 ) xpipe[k] = real(expec);          // standard i/o
      else if( pipe[k] == 4*i+2 ) xpipe[k] = imag(expec);
      else if( pipe[k] == 4*i+3 ) xpipe[k] = real(expec2);
      else if( pipe[k] == 4*i+4 ) xpipe[k] = imag(expec2);
    }
  }
  if( nX > 0 )
    printf("%lG %lG %lG %lG %lG %d %d\n",
      t,xpipe[0],xpipe[1],xpipe[2],xpipe[3],psi.size(),(*stepper).getNumdts());

  for (int n=1; n<=nOfSteps; n++) {
    for (int n1=0; n1<dtsPerStep; n1++) {
      (*stepper)(psi,t,dt,dtlast,rndm);
      t += dt;
      for(int k=0; k<move; k++) {		// use moving basis & cutoff
        psi.adjustCutoff(k,delta,width);
        psi.recenter(k,moveEps);
      }
    }
    for( int i=0; i<nX; i++ ) {		// compute expectation values for
      psi1 = psi;			// output...
      psi1 *= X[i];
      Complex expec = psi*psi1;
      psi1 *= X[i];
      Complex expec2 = psi*psi1 - expec*expec;
      fprintf( fp[i], "%lG %lG %lG %lG %lG\n", 
	      t, real(expec),imag(expec), real(expec2),imag(expec2) ); // file
      for( int k=0; k<4; k++ ) {
	if( pipe[k] == 4*i+1 ) xpipe[k] = real(expec);          // standard i/o
	else if( pipe[k] == 4*i+2 ) xpipe[k] = imag(expec);
	else if( pipe[k] == 4*i+3 ) xpipe[k] = real(expec2);
	else if( pipe[k] == 4*i+4 ) xpipe[k] = imag(expec2);
      }
    }
    if( nX > 0 )
      printf("%lG %lG %lG %lG %lG %d %d\n",
        t,xpipe[0],xpipe[1],xpipe[2],xpipe[3],psi.size(),(*stepper).getNumdts());
  }

  if( nX > 0 ) {
    for ( int i=0; i<nX; i++) fclose(fp[i]); 
    delete fp;
  }
  // Save state to output file
  if (savedState) {
    ofstream os(savedState, ios::out );
    if (!os) {
      cerr << "Can't open output file: " << savedState << endl;
      exit(1);
    }
    os << psi;
  }
}

void Trajectory::sumExp( int nX, const Operator* X, char** fname,
	int dtsPerStep, int nOfSteps, 
        int nTrajectory, int nTrajSave, int ReadFile,
        int move, double delta, int width, double moveEps)
//
// This member function calculates nTrajectory trajectories consisting of 
// nOfSteps output steps, and at each step stores the expectation values 
// of the nX operators X.
// Every nTrajSave trajectory, the mean is computed and saved into files.
// If ReadFile=1 then the files `fname[i]' are read at the beginig of the 
// computation and their values take into account for the mean computation.  
// The variable move contains the number of freedoms for which the moving 
// basis and changing cutoff algorithms should be used;
// delta is the probability threshold for the cutoff adjustment, and
// width is the size of the "pad" for the cutoff.  
//
{
  double t=t0;				// time
  double dtlast=dt;			// keep track of numerical timestep
  State psi0 = psi;                     // initial state
  State psi1 = psi;		        // temporary state
  int nTrajEff = 0;                     // number of trajectories
  FILE** fp;
  fp = new FILE*[nX];                   // nX files are used
  int i,j,k,n;

  // Intialisation of the table OutputExp 
  double** OutputExpR;                  // table of real expectations values
  double** OutputExpI;                  // table of imag expectations values
  OutputExpR = new double*[nOfSteps+1];
  OutputExpI = new double*[nOfSteps+1];
  for ( j=0; j<=nOfSteps; j++) OutputExpR[j] = new double[2*nX];
  for ( j=0; j<=nOfSteps; j++) OutputExpI[j] = new double[2*nX];
  for ( n=0; n<=nOfSteps; n++)
    for ( i=0; i<2*nX; i++) OutputExpR[n][i] = 0.0;
  for ( n=0; n<=nOfSteps; n++)
    for ( i=0; i<2*nX; i++) OutputExpI[n][i] = 0.0;

  //
  // Read the previous results
  //
  if (ReadFile) {
    // open the files in read mode
    for ( i=0; i<nX; i++) fp[i] = fopen(fname[i],"r");
    // check for the files open
    int OpenAllFiles = 1;
    for ( i=0; i<nX; i++) if (fp[i]==NULL) OpenAllFiles=0;
    if (OpenAllFiles) {
      // number of trajectories already computed
      char dummy[100];
      for ( i=0; i<nX; i++) 
	fscanf( fp[i], "%s %d ", dummy, &nTrajEff );
      // read the preceding expectations values
      for ( n=0; n<=nOfSteps; n++)
	for ( i=0; i<nX; i++)
	  fscanf( fp[i], "%lf %lf %lf %lf %lf", 
		 &OutputExpR[n][i], &OutputExpI[n][i], 
		 &OutputExpR[n][nX+i], &OutputExpI[n][nX+i]);
      // close the files
      for ( i=0; i<nX; i++) fclose(fp[i]);
      // transform mean values  
      for ( n=0; n<=nOfSteps; n++)
	for ( i=0; i<2*nX; i++) OutputExpR[n][i] *= nTrajEff;
      for ( n=0; n<=nOfSteps; n++)
	for ( i=0; i<2*nX; i++) OutputExpI[n][i] *= nTrajEff;
      cout << " Using " << nTrajEff << " previous trajectories " << endl;
    }
  }
  //
  // Computation of the expectation values 
  //
  for( k=0; k<nTrajectory; k++ ) {
    cout << "Computing trajectory " << k+1  << endl;
    ++nTrajEff;
    t = t0;
    psi = psi0;

    n=0; 
    for( i=0; i<nX; i++ ) {		// compute expectation values for
      psi1 = psi;			// output...
      psi1 *= X[i];
      Complex expec = psi*psi1;
      psi1 *= X[i];
      Complex expec2 = psi*psi1 - expec*expec;
      OutputExpR[n][i] += real(expec);
      OutputExpI[n][i] += imag(expec);
      OutputExpR[n][nX+i] += real(expec2);
      OutputExpI[n][nX+i] += imag(expec2);
    }

    for ( n=1; n<=nOfSteps; n++) {
      for (int n1=0; n1<dtsPerStep; n1++) {
        (*stepper)(psi,t,dt,dtlast,rndm);
        t += dt;
        for( j=0; j<move; j++) {		// use moving basis & cutoff
          psi.adjustCutoff(j,delta,width);
          psi.recenter(j,moveEps);
        }
      }
      for( i=0; i<nX; i++ ) {		// compute expectation values for
        psi1 = psi;			// output...
        psi1 *= X[i];
        Complex expec = psi*psi1;
        psi1 *= X[i];
        Complex expec2 = psi*psi1 - expec*expec;
        OutputExpR[n][i] += real(expec);
        OutputExpI[n][i] += imag(expec);
        OutputExpR[n][nX+i] += real(expec2);
        OutputExpI[n][nX+i] += imag(expec2);
      }
    }
    //
    // Temporary save of results every nTrajSave 
    //
    if ((nTrajSave>0) && ((k+1) % nTrajSave)==0) {
      for ( i=0; i<nX; i++) fp[i] = fopen(fname[i],"w"); 
      t=t0;
      for ( i=0; i<nX; i++) 
        fprintf( fp[i], "Number_of_Trajectories  %d \n", nTrajEff);
      for ( n=0; n<=nOfSteps; n++) {
        for ( i=0; i<nX; i++) 
          fprintf( fp[i], "%lG %lG %lG %lG %lG\n", t, 
		  OutputExpR[n][i]/nTrajEff, OutputExpI[n][i]/nTrajEff, 
		  OutputExpR[n][nX+i]/nTrajEff, OutputExpI[n][nX+i]/nTrajEff);
	t += dtsPerStep*dt;
      } 
      for ( i=0; i<nX; i++) fclose(fp[i]); 
    }
  }  // End of For Traj Loop

  // Computation of the mean value
  for ( n=0; n<=nOfSteps; n++)
    for ( i=0; i<2*nX; i++) OutputExpR[n][i] /= nTrajEff;
  for ( n=0; n<=nOfSteps; n++)
    for ( i=0; i<2*nX; i++) OutputExpI[n][i] /= nTrajEff;

  //  
  // Writing the results in the files
  //
  for ( i=0; i<nX; i++) fp[i] = fopen(fname[i],"w"); 
  t=t0;
  for ( i=0; i<nX; i++) 
    fprintf( fp[i], "Number_of_Trajectories  %d \n", nTrajEff);
  for ( n=0; n<=nOfSteps; n++) {
    for ( i=0; i<nX; i++) 
      fprintf( fp[i], "%lG %lG %lG %lG %lG\n", t, 
	      OutputExpR[n][i], OutputExpI[n][i], 
	      OutputExpR[n][nX+i], OutputExpI[n][nX+i]);
    t += dtsPerStep*dt;
  } 
  for ( i=0; i<nX; i++) fclose(fp[i]); 
  delete fp;

  // Destruction of the table OutputExp
  for ( j=0; j<nOfSteps; j++) delete OutputExpR[j];
  for ( j=0; j<nOfSteps; j++) delete OutputExpI[j];
  delete OutputExpR;
  delete OutputExpI;
}

Order2Step::Order2Step(const State& psi, const Operator& theH, int theNL,
	const Operator* theL)
{
  H = theH;				// store private data
  nL = theNL;
  L = new Operator[nL];
  Ldag = new Operator[nL];
  dxi = new Complex[nL];
  for (int i=0; i<nL; i++) {
    L[i] = theL[i];
    Ldag[i] = L[i].hc();
    dxi[i] = 0;
  }
  stochasticFlag=0;
  numdtsUsed=0;

  psi2=psi;				// initialize temporaries
  dpsi=psi; 
}

void Order2Step::operator()(State& psi, double t, double dt,
	double& dtlast, ComplexRandom* rndm)
//
// Second order integration of the deterministic part: 
//      see Numerical Recipes, 2nd edition, Eq. (16.1.2).
// Euler for the stochastic part.
{
  Complex lexpect;
 
   
  if( rndm != 0 ) {				// any noise?
    double sqrtdt = sqrt(dt);
    stochasticFlag = 1;
	// set flag for derivs to generate stoch. part and store in newsum

    for (int i=0; i<nL; i++) {			// for each Lindblad...
      dxi[i] = sqrtdt*(*rndm)();		// ...independent noise
    }
  }
  newsum = psi;					// initialize temp
  newsum = 0;
  derivs(t,psi,dpsi);				// also calc. stoch. terms
  dpsi *= (0.5*dt);
  // dpsi = k_1/2 = h/2 f(x_n,y_n) (cf. Numerical Recipes)
  psi2 = psi;
  psi2 += dpsi;       // psi2 = y_n + k_1/2 
  derivs(t+0.5*dt,psi2,dpsi);
  dpsi *= dt;
  psi += dpsi;      // dpsi = k_2 (cf. Numerical Recipes)
  psi += newsum;    // newsum = (Stoch. part)
  psi.normalize();
  numdtsUsed++;
}

Order4Step::Order4Step(const State& psi, const Operator& theH, int theNL,
	const Operator* theL)
{
  H = theH;				// store private data
  nL = theNL;
  L = new Operator[nL];
  Ldag = new Operator[nL];
  dxi = new Complex[nL];
  for (int i=0; i<nL; i++) {
    L[i] = theL[i];
    Ldag[i] = L[i].hc();
    dxi[i] = 0;
  }
  stochasticFlag=0;
  numdtsUsed=0;

  psi0=psi; 				// initialize temporaries
  psi2=psi; 
  dpsi=psi;
}

void Order4Step::operator()(State& psi, double t, double dt,
	double& dtlast, ComplexRandom* rndm)
//
// Fourth order integration of the deterministic part: 
//      see Numerical Recipes, 2nd edition, Eq. (16.1.3).
// Euler for the stochastic part.
// psi is replaced by the new state. psi1 is temporary storage.
{
  Complex lexpect;
  
  if( rndm != 0 ) {				// any noise?
    double sqrtdt = sqrt(dt);
    stochasticFlag = 1;
	// set flag for derivs to generate stoch. part and store in newsum

    for (int i=0; i<nL; i++) {			// for each Lindblad...
      dxi[i] = sqrtdt*(*rndm)();		// ...independent noise
    }
  }
  newsum = psi;					// initialize temp
  newsum = 0;
  psi0 = psi;
  derivs(t,psi,dpsi);			// Also calc. stoch. terms
  dpsi *= (dt/6);
  psi += dpsi;          // dpsi = k1/6

  dpsi *= 3;
  psi2 = psi0;
  psi2 += dpsi;          // psi2 = psi0 + k1/2
  derivs(t+0.5*dt,psi2,dpsi);
  dpsi *= (dt/3);
  psi += dpsi;           // dpsi = k2/3

  dpsi *= 1.5;
  psi2 = psi0;
  psi2 += dpsi;          // psi2 = psi0 + k2/2
  derivs(t+0.5*dt,psi2,dpsi);
  dpsi *= (dt/3);
  psi += dpsi;           // dpsi = k3/3

  dpsi *= 3;
  psi2 = psi0;
  psi2 += dpsi;           // psi2 = psi0 + k3
  derivs(t+0.5*dt,psi2,dpsi);
  dpsi *= (dt/6);
  psi += dpsi;           // psi += k4/6
  psi += newsum;         // psi += (Stoch. part)
  psi.normalize();
  numdtsUsed++;
}

int IntegrationStep::getNumdts()
//
// Return number of deterministic timesteps performed since last
// call to getNumdts
{
  int i;
  i = numdtsUsed;
  numdtsUsed = 0;
  return i;
}

void AdaptiveStep::rkck(double t, double h)
//
// Given values for a state y and its derivative dydt known at t, use the
// fifth-order Cash-Karp Runge-Kutta method to advance the solution over an
// interval h and return the incremented state as yout.  Also return an
// estimate of the local truncation error in yout using the embedded
// fourth-order method.  The user supplies the routine derivs(t,y,dydt) which
// returns the derivative dydt at t.
//
// Adapted from the Numerical Recipes routine RKCK.C
//
{
  static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
	b31=3.0/8.0,b32=9.0/40.0,b41=4.0,b42 = -4.0,b43=1.2,
	b51 = -55.0/81.0, b52=-25.0/9.0, b53 = -175.0/81.0,b54=35.0/27.0,
	b61=-1631.0/11264.0,b62=35.0/256.0,b63=-115.0/7168.0,
	b64=1265.0/4096.0,b65=253.0/4096.0,c1=37888.0/11417.0,
	c3=5120.0/529.0,c4=10240.0/19481.0,c6=512.0/1771.0,
	dc5 = -554.0/1771.0;
  double dc1=-831.0/18944.0,dc3=831.0/17920.0,
	dc4=-831.0/5120.0,dc6=277.0/2048.0;

  ytemp = y;
  derivs(t,ytemp,dydt);
  dydt *= (h*b21);
  ytemp += dydt;				// First step
  derivs(t+a2*h,ytemp,ak2);			// Second step
  ytemp = y;
  dydt *= b31;
  ak2 *= (h*b32);
  ytemp += dydt;
  ytemp += ak2;
  derivs(t+a3*h,ytemp,ak3);			// Third step
  dydt *= b41;
  ak2 *= b42;
  ak3 *= (h*b43);
  ytemp = y;
  ytemp += dydt;
  ytemp += ak2;
  ytemp += ak3;
  derivs(t+a4*h,ytemp,ak4);			// Fourth step
  dydt *= b51;
  ak2 *= b52;
  ak3 *= b53;
  ak4 *= (h*b54);
  ytemp = y;
  ytemp += dydt;
  ytemp += ak2;
  ytemp += ak3;
  ytemp += ak4;
  derivs(t+a5*h,ytemp,ak5);			// Fifth step
  dydt *= b61;
  ak2 *= b62;
  ak3 *= b63;
  ak4 *= b64;
  ak5 *= (h*b65);
  ytemp = y;
  ytemp += dydt;
  ytemp += ak2;
  ytemp += ak3;
  ytemp += ak4;
  ytemp += ak5;
  derivs(t+a6*h,ytemp,ak2);			// Sixth step
//
// Accumulate increments with proper weights
//
  dydt *= c1;
  ak3 *= c3;
  ak4 *= c4;
  ak2 *= (h*c6);
  yout = y;
  yout += dydt;
  yout += ak3;
  yout += ak4;
  yout += ak2;
//
// Estimate error as difference between fourth and fifth order methods.
//
  dydt *= dc1;
  ak3 *= dc3;
  ak4 *= dc4;
  ak5 *= dc5;
  ak2 *= dc6;
  yerr = dydt;
  yerr += ak3;
  yerr += ak4;
  yerr += ak5;
  yerr += ak2;
}

void AdaptiveStep::rkqs(double& t, double htry, double eps,
	double& hdid, double& hnext)
//
// Fifth-order Runge-Kutta step with monitoring of local truncation error to
// ensure accuracy and adjust stepsize.  Input are the State y and its
// derivative dydt at the starting value of the time t.  Also input are the
// stepsize to be attempted, htry, and the required accuracy, eps. On output,
// y and t are replaced by their new values, hdid is the stepsize that was
// actually accomplished, and hnext is the estimated next stepsize.  derivs
// is the user-supplied routine that computes the right-hand side derivative.
//
// Adapted from the Numerical Recipes routine RKQS.C
//
{
  double errmax,h,tnew;

  h=htry;
  for (;;) {
	rkck(t,h);
	errmax = sqrt(real(yerr*yerr));
	errmax /= eps;
	if (errmax > 1.0) {
		h=SAFETY*h*pow(errmax,PSHRNK);
//		if (h < 0.1*h) h *= 0.1;
		tnew=t+h;
		if (tnew == t) {
		  cerr << "stepsize underflow in rkqs" << endl;
		  cerr << "errmax = " << errmax << endl;
		  cerr << "h = " << h << endl;
		  cerr << "eps = " << eps << endl;
		  exit(1);
		}
		continue;
	} else {
		if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
		else hnext=5.0*h;
		t += (hdid=h);
		y = yout;
		numdtsUsed++;
		break;
	}
  }
}


void AdaptiveStep::odeint(State& ystart, double t1, double t2, double eps, double& h1,
	double hmin, int& nok, int& nbad, Complex* dxi)
//
// Runge-Kutta driver with adaptive stepsize control  Integrate starting values
// ystart from t1 to t2 with accuracy eps.  h1 should be set as a guessed
// first stepsize, hmin as the minimum allowed stepsize (can be zero).  On
// output nok and nbad are the number of good and bad (but retried and fixed)
// steps taken, and ystart is replaced by values at the end of the integration
// interval.  derivs is the routine for calculating the right-hand side of
// the QSD equation, while rkqs is the name of the stepper routine
// to be used.
//
// Adapted from Numerical Recipes routine ODEINT.C
//
{
  int nstp;
  double t,hnext,hdid,h;

  t = t1;
  h = h1; if( h > (t2-t1) ) h = (t2-t1);
  nok = nbad = 0;
  y = ystart;
  for (nstp=1;nstp<=MAXSTP;nstp++) {
	if ((t+h-t2)*(t+h-t1) > 0.0) h=t2-t;
	rkqs(t,h,eps,hdid,hnext);		// deterministic part
	if (hdid == h) ++nok; else ++nbad;
	y.normalize();				// renormalize
	if ((t-t2)*(t2-t1) >= 0.0) {
		ystart = y;
		h1 = hnext;
		if( listPtr < listDim ) {	// if keeping track...
		  dtNumList[listPtr] = nok+nbad;  // ...store number of steps
		  listPtr++;
		}
		return;
	}
	if (fabs(hnext) <= hmin) {
	  cerr << "Step size too small in odeint" << endl;
	  exit(1);
	}
	h=hnext;
  }
  cerr << "Too many steps in routine odeint" << endl;
  exit(1);
}

void AdaptiveStochStep::odeint(State& ystart, double t1, double t2, double eps, double& h1,
	double hmin, int& nok, int& nbad, ComplexRandom* rndm) 
//
// Runge-Kutta driver with adaptive stepsize control  Integrate starting values
// ystart from t1 to t2 with accuracy eps.  h1 should be set as a guessed
// first stepsize, hmin as the minimum allowed stepsize (can be zero).  On
// output nok and nbad are the number of good and bad (but retried and fixed)
// steps taken, and ystart is replaced by values at the end of the integration
// interval.  derivs is the routine for calculating the right-hand side of
// the QSD equation, while rkqs is the name of the stepper routine
// to be used.
//
// Also, generate one stochastic timestep per deterministic timestep;
//
// Adapted from Numerical Recipes routine ODEINT.C
//
{
  int nstp;
  double t,hnext,hdid,h;
  double sqrtdt;

  t = t1;
  h = h1; if( h > (t2-t1) ) h = (t2-t1);
  nok = nbad = 0;
  y = ystart;
  for (nstp=1;nstp<=MAXSTP;nstp++) {
    if ((t+h-t2)*(t+h-t1) > 0.0) h=t2-t;
    if( rndm != 0 ) {				// any noise?
      stochasticFlag = 1;
	  // set flag for derivs to generate stoch. part and store in newsum
      for (int i=0; i<nL; i++)			// for each Lindblad...
        dxi[i] = (*rndm)();			// ...independent noise
    }
    newsum = 0;
    rkqs(t,h,eps,hdid,hnext);			// deterministic part
    if (hdid == h) ++nok; else ++nbad;
    if( rndm != 0 ) {
      sqrtdt = sqrt(hdid);
      newsum *= sqrtdt;
      y += newsum;
    }
    y.normalize();				// renormalize
    if ((t-t2)*(t2-t1) >= 0.0) {
      ystart = y;
      h1 = hnext;
      if( listPtr < listDim ) {			// if keeping track...
        dtNumList[listPtr] = nok+nbad;		// ...store number of steps
        listPtr++;
      }
      return;
    }
    if (fabs(hnext) <= hmin) {
      cerr << "Step size too small in odeint" << endl;
      exit(1);
    }
    h=hnext;
  }
  cerr << "Too many steps in routine odeint" << endl;
  exit(1);
}


void AdaptiveJump::odeint(State& ystart, double t1, double t2, double eps, double& h1,
	double hmin, int& nok, int& nbad, ComplexRandom* rndm) 
//
// Runge-Kutta driver with adaptive stepsize control  Integrate starting values
// ystart from t1 to t2 with accuracy eps.  h1 should be set as a guessed
// first stepsize, hmin as the minimum allowed stepsize (can be zero).  On
// output nok and nbad are the number of good and bad (but retried and fixed)
// steps taken, and ystart is replaced by values at the end of the integration
// interval.  derivs is the routine for calculating the right-hand side of
// the QSD equation, while rkqs is the name of the stepper routine
// to be used.
//
// Also, generate one stochastic timestep per deterministic timestep;
//
// Adapted from Numerical Recipes routine ODEINT.C
//
{
  int nstp,count;
  double t,hnext,hdid,h;
  double prob,prob2,total;

  t = t1;
  h = h1; if( h > (t2-t1) ) h = (t2-t1);
  nok = nbad = 0;
  y = ystart;
  for (nstp=1;nstp<=MAXSTP;nstp++) {
    if ((t+h-t2)*(t+h-t1) > 0.0) h=t2-t;
    if( rndm != 0 ) {				// any noise?
      stochasticFlag = 1;			// have derivs calc <Ldag L>
    }
    newsum = y;                                 // keep the actual value of psi
    rkqs(t,h,eps,hdid,hnext);			// deterministic part
    if (hdid == h) ++nok; else ++nbad;
    if( rndm != 0 ) {				// any noise?
      total = 0;
      for(count=0; count<nL; count++)
	total += lindVal[count];		// sum up <Ldag L>
// Two different choice for the test-jump
// they should give the same results
      prob = 1.0 - real(y*y);
//      prob = total*hdid;
//
      if( real((*rndm)()) < prob ) {		// make a jump?
	y = newsum;
        prob = real((*rndm)());
	for(count=0; prob>0; count++) {
          prob2 = lindVal[count]/total;
          if( prob < prob2 ) {			// jump !
            y *= L[count];			// L|psy>
            if( outFile != 0 )			// output jump time to file
              fprintf(outFile,"%lf  %d\n",t,count);
          }
          prob -= prob2;
        }
	if (t2-t > hnext) hnext = t2-t;  // After a jump the step h
	                                 // is re-initialized 
      }
    }
    y.normalize();				// renormalize
    if ((t-t2)*(t2-t1) >= 0.0) {
      ystart = y;
      h1 = hnext;
      if( listPtr < listDim ) {			// if keeping track...
        dtNumList[listPtr] = nok+nbad;		// ...store number of steps
        listPtr++;
      }
      return;
    }
    if (fabs(hnext) <= hmin) {
      cerr << "Step size too small in odeint" << endl;
      exit(1);
    }
    h=hnext;
  }
  cerr << "Too many steps in routine odeint" << endl;
  exit(1);
}

void AdaptiveOrthoJump::odeint(State& ystart, double t1, double t2, double eps, double& h1,
	double hmin, int& nok, int& nbad, ComplexRandom* rndm) 
//
// Runge-Kutta driver with adaptive stepsize control  Integrate starting values
// ystart from t1 to t2 with accuracy eps.  h1 should be set as a guessed
// first stepsize, hmin as the minimum allowed stepsize (can be zero).  On
// output nok and nbad are the number of good and bad (but retried and fixed)
// steps taken, and ystart is replaced by values at the end of the integration
// interval.  derivs is the routine for calculating the right-hand side of
// the orthogonal jump equation, while rkqs is the name of the stepper routine
// to be used.
//
// Also, generate one stochastic timestep per deterministic timestep;
//
// Adapted from Numerical Recipes routine ODEINT.C
//
{
  int nstp,count;
  double t,hnext,hdid,h;
  double prob,prob2,total;
  Complex expect;

  t = t1;
  h = h1; if( h > (t2-t1) ) h = (t2-t1);
  nok = nbad = 0;
  y = ystart;
  for (nstp=1;nstp<=MAXSTP;nstp++) {
    if ((t+h-t2)*(t+h-t1) > 0.0) h=t2-t;
    if( rndm != 0 ) {				// any noise?
      stochasticFlag = 1;			// have derivs calc <Ldag L>
    }
    newsum = y;
    rkqs(t,h,eps,hdid,hnext);			// deterministic part
    if (hdid == h) ++nok; else ++nbad;
    if( rndm != 0 ) {				// any noise?
      total = 0;
      for(count=0; count<nL; count++)
	total += lindVal[count];		// sum up <Ldag L>-<Ldag><L>
      prob = total*hdid;
      if( real((*rndm)()) < prob ) {		// make a jump?
	y = newsum;
        prob = real((*rndm)());
        for(count=0; prob>0; count++) {
          prob2 = lindVal[count]/total;
          if( prob < prob2 ) {		        // jump!
            y *= L[count];	                // L |psi>
	    expect = newsum*y;                  // <L>
	    newsum *= expect;                   // <L>|psi>
	    y -= newsum;                        // L|psi> - <L>|psi>
            if( outFile != 0 )			// output jump time to file
              fprintf(outFile,"%lf  %d\n",t,count);
	  }
          prob -= prob2;
        }
	if (t2-t > hnext) hnext = t2-t;  // After a jump the step h
	                                 // is re-initialized 
      }
    }
    y.normalize();				// renormalize
    if ((t-t2)*(t2-t1) >= 0.0) {
      ystart = y;
      h1 = hnext;
      if( listPtr < listDim ) {			// if keeping track...
        dtNumList[listPtr] = nok+nbad;		// ...store number of steps
        listPtr++;
      }
      return;
    }
    if (fabs(hnext) <= hmin) {
      cerr << "Step size too small in odeint" << endl;
      exit(1);
    }
    h=hnext;
  }
  cerr << "Too many steps in routine odeint" << endl;
  exit(1);
}

void AdaptiveStep::dtListSet( int theDim )
//
// Sets up a list to store the number of deterministic steps per
// stochastic timestep
//
{
#ifndef OPTIMIZE_QSD
  if( theDim < 1 )
    error("Negative dimension specified in AdaptiveStep::dtListSet!");
#endif
  listDim = theDim;
  dtNumList = new int[theDim];
}

void AdaptiveStep::dtListRead()
//
// Prints out the number of deterministic steps per stochastic timestep
//
{
  for( int i=0; i<listPtr; i++ )
    cout << dtNumList[i] << endl;
}

int AdaptiveStep::dtListElem( int theElem )
//
// Prints out the number of deterministic steps in a particular
// stochastic timestep
//
{
#ifndef OPTIMIZE_QSD
  if( (theElem < 0) || (theElem >= listPtr) )
    error("Illegal list element in AdaptiveStep::dtListElem!");
#endif
  return dtNumList[theElem];
}

void AdaptiveStep::dtListClear()
//
// Clear the list of numbers of deterministic steps per stochastic timestep
//
{
  if( listDim != 0 )
#ifndef NON_GNU_DELETE
    delete[] dtNumList;			// deallocate memory
#else
    delete[listDim] dtNumList;		// deallocate memory
#endif
  listDim = 0;
  listPtr = 0;
  dtNumList = 0;
}

void AdaptiveStep::dtListReset()
//
// Reset the array of numbers without deallocating it
//
{
  listPtr = 0;
}

void IntegrationStep::derivs(double t, State& psi, State& dpsi)
//
// Computes the change of the state psi in a timestep dt due to the
// deterministic part of the QSD equation.
//
// Also, when stochasticFlag is set, computes the sum for the stochastic part.
//
{
  Complex expect;
  double sum=0;
  Complex stochSum=0;

  dpsi = psi;
  dpsi *= H(t);
  dpsi *= M_IM;
  for (int count=0; count<nL; count++) {
    temp1 = psi;
    temp1 *= L[count];
    expect = psi*temp1;
    if( stochasticFlag != 0 ) {
      temp0 = temp1;
      temp0 *= dxi[count];
      newsum += temp0;
      stochSum -= expect*dxi[count] + 0.5*conj(expect)*expect*(dxi[count]*dxi[count] - 0.0);
    }
    sum -= norm(expect);
    temp0 = temp1;
    temp0 *= conj(expect);
    dpsi += temp0;
    temp1 *= Ldag[count];
    temp1 *= -0.5;
    dpsi += temp1;
  }
  sum *= 0.5;
  temp0 = psi;
  temp0 *= sum;
  dpsi += temp0;
  if( stochasticFlag != 0 ) {
    temp0 = psi;
    temp0 *= stochSum;
    newsum += temp0;
  }
  stochasticFlag = 0;
}

void AdaptiveJump::derivs(double t, State& psi, State& dpsi)
//
// Computes the change of the state psi in a timestep dt due to the
// deterministic part of the quantum jumps evolution.
//
// Also, when stochasticFlag is set, computes the sum for the stochastic part.
//
{
  dpsi = psi;
  dpsi *= H(t);
  dpsi *= M_IM;                                     // -iH|psi>
  for (int count=0; count<nL; count++) {
    temp1 = psi;
    temp1 *= L[count];                              // L|psi>
    if (stochasticFlag == 1)
      lindVal[count] = real(temp1*temp1);           // <Ldag L>
    temp1 *= Ldag[count];                           // Ldag L|psi>
    temp1 *= -0.5;
    dpsi += temp1;                                  // -0.5 Ldag L|psi>  
  }
  stochasticFlag = 0;
}

void AdaptiveOrthoJump::derivs(double t, State& psi, State& dpsi)
//
// Computes the change of the state psi in a timestep dt due to the
// deterministic part of the orthogonal jump evolution.
//
// Also, when stochasticFlag is set, computes the sum for the stochastic part.
//
{
  Complex expect;

  dpsi = psi;
  dpsi *= H(t);
  dpsi *= M_IM;                 // dpsi = -iH|psi>
  for (int count=0; count<nL; count++) {
    temp1 = psi;
    temp1 *= L[count];          // L|psi>
    expect = psi*temp1;         // <L>=<psi|L|psi>
    if (stochasticFlag == 1)
      lindVal[count] = real(temp1*temp1)-norm(expect);  // <LdagL>-<Ldag><L>
    temp0 = temp1;              
    temp0 *= conj(expect);      // <Ldag>L|psi>
    dpsi += temp0;              // dpsi += <Ldag> L |psi>
    temp1 *= Ldag[count];       // Ldag L |psi>
    temp1 *= -0.5;
    dpsi += temp1;              // dpsi += -0.5 Ldag L |psi> 
  }
  stochasticFlag = 0;
}

void Trajectory::error(char* message)	// print error message and exit
{
  cerr << message << endl;
  exit(1);
}

void IntegrationStep::error(char* message)	// print error message
{						// and exit
  cerr << message << endl;
  exit(1);
}

#undef MAXSTP
#undef TINY
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON




