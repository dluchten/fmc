#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_
#include "grid2d.h"
#include "velocityfield.h"

class Integrator {
public:
  Integrator(const double dt, const VelocityField &vf) 
      : dt_(dt), dimen_(vf.dimen()), vf_(vf) {}
  virtual ~Integrator() {}
  virtual void Step(const double &t, double *x) = 0;
  void Step(const double &t, const Grid2D &gridin, Grid2D &gridout);
protected:
  const double dt_;                     // timestep
  const int dimen_;                     // dimension of state x
  const VelocityField &vf_;             // functor to evaluate f(x,t)
};

#endif  // INTEGRATOR_H_
