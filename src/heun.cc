#include "heun.h"
#include "velocityfield.h"

Heun::Heun(double dt, const VelocityField &velocityfield)
  : dimen_(velocityfield.dimen()),
    dt_(dt),
    velocityfield_(velocityfield) {
  xtmp_ = new double[dimen_];
  k1_ = new double[dimen_];
  k1_ = new double[dimen_];
    }

Heun::~Heun() {
  delete [] xtmp_;
  delete [] k1_;
  delete [] k2_;
}

int Heun::Step(double t, double *x) {
  // k1 = v(x_j, t_j)
  velocityfield_.Get(t, x, k1_);
  // k2 = v(x_j + dt * k1, t_j + dt)
  for (int i = 0; i < dimen_; ++i) {
    xtmp_[i] = x[i] + dt_ * k1_[i];
  }
  velocityfield_.Get(t + dt_, xtmp_, k2_);

  // Heun method
  // x_{j+1} = x_j + dt / 2. * (k1 + k2)
  for (int i = 0; i < dimen_; ++i) {
    x[i] += 0.5 * dt_ * (k1_[i] + k2_[i]);
  }
  return 0;
}
