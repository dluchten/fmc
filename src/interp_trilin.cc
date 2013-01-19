#include "interp_trilin.h"

InterpTrilin::InterpTrilin(const VecDoub1D &x0, const VecDoub1D &x1, const VecDoub1D &x2, VecDoub3D *fx)
  : Interp(3),
    x0_(x0),
    x1_(x1),
    x2_(x2),
    fx_(fx) {
  stencil_.push_back(new Stencil_(x0, 2));
  stencil_.push_back(new Stencil_(x1, 2));
  stencil_.push_back(new Stencil_(x2, 2));
  low_edge_ = new int[3];
    }

InterpTrilin::~InterpTrilin() {
  while (!stencil_.empty()) {
    delete stencil_.back();
    stencil_.pop_back();
  }
  delete [] low_edge_;
}

double InterpTrilin::RawInterpolate_(const double *x) const {
  int i = low_edge_[0];
  int j = low_edge_[1];
  int k = low_edge_[2];

  double m0, m1, m2, fx;

  // Weigh by position along each axis
  m0 = (x[0]-x0_[i]) / (x0_[i+1]-x0_[i]);
  m1 = (x[1]-x1_[j]) / (x1_[j+1]-x1_[j]);
  m2 = (x[2]-x2_[k]) / (x2_[k+1]-x2_[k]);

  // Interpolate
  fx = (1.-m0) * (1.-m1) * (1.-m2) * (*fx_)[i  ][j  ][k  ]
     + (1.-m0) * (1.-m1) *     m2  * (*fx_)[i  ][j  ][k+1]
     + (1.-m0) *     m1  * (1.-m2) * (*fx_)[i  ][j+1][k  ]
     + (1.-m0) *     m1  *     m2  * (*fx_)[i  ][j+1][k+1]
     +     m0  * (1.-m1) * (1.-m2) * (*fx_)[i+1][j  ][k  ]
     +     m0  * (1.-m1) *     m2  * (*fx_)[i+1][j  ][k+1]
     +     m0  *     m1  * (1.-m2) * (*fx_)[i+1][j+1][k  ]
     +     m0  *     m1  *     m2  * (*fx_)[i+1][j+1][k+1];

  return fx;
}
