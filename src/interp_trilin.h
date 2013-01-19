#ifndef INTERP_TRILIN_H_
#define INTERP_TRILIN_H_

#include "interp.h"

class InterpTrilin : public Interp {
 public:
  InterpTrilin(const VecDoub1D &x0, const VecDoub1D &x1, const VecDoub1D &x2, VecDoub3D *fx);
  ~InterpTrilin();
  void ChangeFieldPtr(VecDoub3D *ptr) { fx_ = ptr; }
 private:
  double RawInterpolate_(const double *x) const;
  const VecDoub1D &x0_;
  const VecDoub1D &x1_;
  const VecDoub1D &x2_;
  VecDoub3D *fx_;
};

#endif  // INTERP_TRILIN_H_
