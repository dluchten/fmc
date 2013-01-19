#include "integrator.h"

void Integrator::Step(const double &t, const Grid2D &gridin, Grid2D &gridout) {
  assert(dimen_ == 2);
  double x[2];
  int m = (int)gridin.Length1();
  int n = (int)gridin.Length2();
  for (int i = 0; i < m; ++i) {
    for (int j =  0; j < n; ++j) {
      // Get (x, y)
      x[0] = gridin.X(i,j); // x coord
      x[1] = gridin.Y(i,j); // y coord
      // integrate for one time step
      Step(t, x);
      // Set new (x, y)
      gridout.X(i,j) = x[0];
      gridout.Y(i,j) = x[1];
    }
  }
}      
