#include "doublegyre.h"
#include <cmath>
using namespace std;

// Static members need to be pre-declared
const double DoubleGyre::kPi = acos(-1.);

DoubleGyre::DoubleGyre(double A, double eps, double omega)
  : A_(A = 0.1),
    eps_(eps = 0.25), 
    omega_(omega = 2. * kPi)
{}

DoubleGyre::~DoubleGyre() {}

int DoubleGyre::velocity_field(double t, const double *x, double *vf) const {
  double a = eps_ * sin(omega_ * t);
  double b = 1 - 2. * a;
  double f = a * x[0] * x[0] + b * x[0];
  double df = 2. * a * x[0] + b;

  vf[0] = -kPi * A_ * sin(kPi * f) * cos(kPi * x[1]);
  vf[1] =  kPi * A_ * cos(kPi * f) * sin(kPi * x[1]) * df;

  return 0;
}
