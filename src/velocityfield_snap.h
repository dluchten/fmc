#ifndef VELOCITYFIELD_SNAP_H_
#define VELOCITYFIELD_SNAP_H_

#include "typedefs.h"
#include "velocityfield.h"

class VelocityFieldSnap : public VelocityField {
 public:
 VelocityFieldSnap(const int dimen) : kDimen_(dimen) {};
  virtual ~VelocityFieldSnap() {};
  virtual void Get(double t, const double *x, double *vf) const = 0;
  virtual int Set(const int t, const int *x,
		  const int comp, const double vcomp) = 0;
  int dimen() const { return kDimen_; }
  int num_snaps() const { return t_->size(); }
  double tlo() const { return (*t_)[0]; }
  double thi() const { return t_->back(); }
 protected:
  VecDoub1D *t_;      // vector of time values
  const int kDimen_;  // spatial dimensionality
};

#endif  // VELOCITYFIELD_SNAP_H_
