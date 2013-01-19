#ifndef INTERP_H_
#define INTERP_H_

#include "typedefs.h"

// Base class for multi-dimensional interpolators
// Possible to change the field over which it interpolates
// (useful for vector fields)
// The implementation here follows that in Numerical Recipes in C++, 3rd ed

class Interp {

 public:
  Interp(const int dimen);
  virtual ~Interp() {}
  double Interpolate(const double *x) {
    for (int i=0; i<kDimen_; ++i) {
      low_edge_[i] = (*(stencil_[i])).low_edge(x[i]);
    }
    return RawInterpolate_(x);
  }

 protected:

  // A class that finds the grid points near a given coordinate
  // These grid points are then used in RawInterpolate_
  class Stencil_ {
  public:
    ~Stencil_();
    Stencil_(const VecDoub1D &x, const int size);

    // Return low-index edge of the series of points used for interpolation
    int low_edge(const double x) {
      return hunt_or_loc_ ? Hunt_(x) : Locate_(x);
    }
  private:
    const VecDoub1D &x_;         // array of ordered grid points
    const int kStencilSize_;     // number of elements the stencil looks at
    const int kBias_;            // used in selection of Hunt_ or Locate_
    int hunt_or_loc_;            // selects whether to use Hunt_ or Locate_
    int last_pos_;               // position of stencil from last search
    int Hunt_(const double x);   // search in vicinity of previously found pt
    int Locate_(const double x); // use repeated bisection to find grid pt
  };

  const int kDimen_;               // spatial dimensionality
  std::vector<Stencil_*> stencil_; // one Stencil_ for each dimension
  int *low_edge_;                  // store low-index edge of the stencil

  // Uses input double value and indices found by Stencil_
  // to interpolate over a small number of grid points
  // Different interpolators (linear, cubic) do this differently
  virtual double RawInterpolate_(const double *x) const = 0;
};

#endif  // INTERP_H_
