// Heun method

#ifndef HEUN_H_
#define HEUN_H_

#include "integrator.h"
class VelocityField;

class Heun : public Integrator {
public:
  Heun(double dt, const VelocityField &velocityfield);
  ~Heun();
  int Step(double t, double *x);
private:
  const int dimen_;                     // dimension of state x
  const double dt_;                     // timestep
  const VelocityField &velocityfield_; 
  double *xtmp_;                          
  double *k1_;
  double *k2_;
};

#endif  // HEUN_H_
