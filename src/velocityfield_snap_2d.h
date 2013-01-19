#ifndef VELOCITYFIELD_SNAP_2D_H_
#define VELOCITYFIELD_SNAP_2D_H_

#include "velocityfield_snap.h"
#include "interp_trilin.h"
#include <string>

class VelocityFieldSnap2D : public VelocityFieldSnap {
 public:
  VelocityFieldSnap2D(const VecDoub1D &t,
		      const VecDoub1D &x, const VecDoub1D &y);
  VelocityFieldSnap2D(const VecDoub1D &t, const std::vector<std::string> &flist);
  ~VelocityFieldSnap2D();
  void Get(double t, const double *x, double *vf) const;
  int Set(const int t, const int *x, const int comp, const double vcomp);
  int nx() const { return x_->size(); }
  int ny() const { return y_->size(); }
  double xlo() const { return (*x_)[0]; }
  double ylo() const { return (*y_)[0]; }
  double xhi() const { return x_->back(); }
  double yhi() const { return y_->back(); }
 private:
  VecDoub1D *x_, *y_;        // arrays of spatial coords for 2D
  VecDoub3D *velx_, *vely_;  // gridded x,y velocities
  InterpTrilin *interp_;     // trilinear interpolator
  double *coords_;           // temp space for coord array used by interpolator
};

#endif  // VELOCITYFIELD_SNAP_2D_H_
