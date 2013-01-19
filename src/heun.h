// Heun method

#ifndef HEUN_H_
#define HEUN_H_

#include "integrator.h"

class Heun : public Integrator {
public:
  Heun(const double dt, const VelocityField &vf);
  ~Heun();
  void Step(const double &t, double *x);
private:
  double *xtmp_;                          
  double *k1_;
  double *k2_;
};

#endif  // HEUN_H_
