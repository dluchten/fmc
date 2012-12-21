#ifndef MODEL_H_
#define MODEL_H_

class Model {
public:
  virtual ~Model() {} // MUST have a function body

  // vf = v(x,t)
  virtual int velocity_field(double t, const double *x, double *vf) const = 0;

  // number of states (size of x)
  virtual int dimen() const = 0;
};

#endif  //  MODEL_H_
