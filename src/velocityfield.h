#ifndef VELOCITYFIELD_H_
#define VELOCITYFIELD_H_

class VelocityField {
public:
  virtual ~VelocityField() {}

  // vf = v(x,t)
  virtual void Get(double t, const double *x, double *vf) const = 0;

  // number of states (size of x)
  virtual int dimen() const = 0;
};

#endif  // VELOCITYFIELD_H_
