class Poisson {
// The class Poisson solves the Poisson equation

public:

// constructors and destructors

  Poisson();
  // Defualt constructor; .
  ~Poisson();         // destructor

// public functions used by constructors and destructors

  void poisson_init();
  double* poisson_sum(double* rho);
  double* poisson_fft(double* rho);
  double yukawa_fs(double x);
  int nint(double x);
  double besselint(double x);

  struct fftw
{

  int n1, n2;
  int isReal;        // is the fft real or complex
  fftw_plan planf;        // the plan for forward transform
  fftw_plan planb;        // the plan for backward transform

};

  struct DCF
{
   int n1, n2;
   double* realSpace;
   fftw_complex* fourierSpace;
   int nx; //  = n1/2 + 1, first dimension of the FS array
   fftw* fftp;
};

  DCF cf;

};
