#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

class Integrator {
public:
  virtual ~Integrator() {}
  virtual int Step(const double &t, double *x) = 0;
};

#endif  // INTEGRATOR_H_
