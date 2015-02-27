#include "Complex.h"

class Coeff {

public:

  Coeff();
//  ~Coeff();


void coefficients(double *c);
int  getN();

void coefficients(int n, double *c);

private:

 int m; // derivative order
 double del;
 int n; // order of the approximation
 double* x; // Store for

};
