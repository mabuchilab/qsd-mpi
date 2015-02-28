class Mesh {
// The class Mesh represents the grid on which the calculations are carried out

public:

// constructors and destructors

  Mesh();
  // Defualt constructor; produces a Mesh of size 0.
  Mesh(int n);
  // Produces a Mesh with n by n grid points
  ~Mesh();         // destructor

// public functions used by constructors and destructors

  void laplacian(double*  f, double* lapl, double* c, int i);  // computes the laplacian
 
  int getN();
  double getDelta();

  double dotproduct(double *, double *);
  double x(int ix, int iy);
  double y(int ix, int iy);

};
