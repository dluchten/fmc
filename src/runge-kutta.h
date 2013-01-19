// Fourth-order Runge-Kutta method

#ifndef RUNGEKUTTA_H_
#define RUNGEKUTTA_H_

#include "integrator.h"

class RungeKutta4 : public Integrator {
public:
  RungeKutta4(const double dt, const VelocityField &vf);
  ~RungeKutta4();
  void Step(const double &t, double *x);
private:
  double *xtmp_; // temporary space for R-K integrator
  double *k1_;   // temporary space for R-K integrator
  double *k2_;   // temporary space for R-K integrator
  double *k3_;   // temporary space for R-K integrator
  double *k4_;   // temporary space for R-K integrator
};

#endif  // RUNGEKUTTA_H_
