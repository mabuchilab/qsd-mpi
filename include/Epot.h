class Epot {
// The class Epot represents the external potentials supported

public:

// constructors and destructors

  Epot();
  // Defualt constructor; .
  ~Epot();         // destructor

// public functions used by constructors and destructors

  double* harmonic();
  double* quartic();

};
