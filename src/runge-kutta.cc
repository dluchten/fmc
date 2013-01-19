#include "runge-kutta.h"

RungeKutta4::RungeKutta4(const double dt, const VelocityField &vf)
    : Integrator(dt, vf) {
  xtmp_ = new double[dimen_];
  k1_ = new double[dimen_];
  k2_ = new double[dimen_];
  k3_ = new double[dimen_];
  k4_ = new double[dimen_];
}

RungeKutta4::~RungeKutta4() {
  delete [] xtmp_;
  delete [] k1_;
  delete [] k2_;
  delete [] k3_;
  delete [] k4_;
}

void RungeKutta4::Step(const double &t, double *x) {
  // k1 = f(x_j, t_j)
  vf_.Get(t, x, k1_);
  // k2 = f(x_j + 0.5 * dt * k1, t_j + 0.5 * dt)
  for (int i = 0; i < dimen_; ++i) {
    xtmp_[i] = x[i] + 0.5 * dt_ * k1_[i];
  }
  vf_.Get(t + 0.5 * dt_, xtmp_, k2_);
  // k3 = f(x_j + 0.5 * dt * k2, t_j + 0.5 * dt)
  for (int i = 0; i < dimen_; ++i) {
    xtmp_[i] = x[i] + 0.5 * dt_ * k2_[i];
  }
  vf_.Get(t + 0.5 * dt_, xtmp_, k3_);
  // k4 = f(x_j + dt * k3, t_j + dt)
  for (int i = 0; i < dimen_; ++i) {
    xtmp_[i] = x[i] + dt_ * k3_[i];
  }
  vf_.Get(t + dt_, xtmp_, k4_);

  // Runge Kutta method
  // x_{j+1} = x_j + dt / 6. * (k1 + 2. * (k2 + k3) + k4)
  for (int i = 0; i < dimen_; ++i) {
    x[i] += dt_ / 6. * (k1_[i] + 2. * (k2_[i] + k3_[i]) + k4_[i]);
  }
}
