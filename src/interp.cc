#include "interp.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>

Interp::Interp(const int dimen)
  : kDimen_(dimen) {}

Interp::Stencil_::Stencil_(const VecDoub1D &x, const int size)
  : x_(x),
    kStencilSize_(size),
    kBias_(std::min(1, int(pow(double(x.size()), 0.25)))),
    hunt_or_loc_(0),
    last_pos_(0) {

  int nx = x.size();

  // Check if sizes of stencil cand array are appropriate
  if (nx<2 || kStencilSize_<2 || kStencilSize_>2) {
    std::cout << "Bad size for independent var array or stencil" << std::endl;
    exit(1);
  }

  // Check if x is monotonic
  int first_sign = (x[1] > x[0]) - (x[1] < x[0]);
  int current_sign;
  if (first_sign == 0) {
    std::cout << "Array not monotonic" << std::endl;
    exit(1);
  }
  for (int i=2; i<nx; ++i) {
    current_sign = (x[i] > x[i-1]) - (x[i] < x[i-1]);
    if (current_sign != first_sign) {
      std::cout << "Array not monotonic" << std::endl;
      exit(1);
    }
  }
}

Interp::Stencil_::~Stencil_() {}

int Interp::Stencil_::Hunt_(const double x) {

  int jl = last_pos_;
  int jm, ju;
  int inc = 1;
  int nx = x_.size();
  
  // True if ascending order, false if descending, assuming monotonic
  bool ascnd = (x_[nx-1] >= x_[0]);

  if (jl<0 || jl>nx-1) { // input guess not useful, use bisection
    jl = 0;
    ju = nx-1;
  } else {
    if ((x >= x_[jl]) == ascnd) { // use hunt method
      while(true) {
	ju = jl + inc;
	if (ju >= nx-1) {ju = nx-1; break;}
	else if ((x < x_[ju]) == ascnd) break;
	else {
	  jl = ju;
	  inc += inc;
	}
      }
    } else {
      ju = jl;
      while(true) {
	jl = jl - inc;
        if (jl <= 0) {jl = 0; break;}
        else if ((x >= x_[jl]) == ascnd) break;
        else {
          ju = jl;
          inc += inc;
        }
      }
    }
  }
  while (ju-jl > 1) {
    jm = (ju+jl) >> 1;
    if ((x >= x_[jm]) == ascnd)
      jl = jm;
    else
      ju = jm;
  }
  hunt_or_loc_ = abs(jl-last_pos_) > kBias_ ? 0 : 1;
  last_pos_ = jl;
  return std::max(0, std::min(nx-kStencilSize_, jl-((kStencilSize_-2)>>1)));
}

int Interp::Stencil_::Locate_(const double x) {
  int ju, jm, jl;

  int nx = x_.size();

  // Assume x values ordered monotonically
  bool ascnd = (x_[nx-1] >= x_[0]);

  // Lower and upper cursors start out at [0] and [n-1] of array
  // While the distance between their indices is greater than 1 index,
  //   get the x value of the index between them and compare with given x
  // If x value is >= than mid-index value, and array ascends,
  //   change the low cursor to the current middle index,
  //   otherwise, change the upper cursor to the current middle index.
  jl = 0;
  ju = nx-1;
  while (ju-jl > 1) {
    jm = (ju+jl) >> 1; // take average ... add then shift right (= divide by 2)
    if ((x >= x_[jm]) == ascnd)
      jl = jm;
    else
      ju = jm;
  }

  // If the index distance from the last index starting point (hunt or locate),
  //   is bigger than some value, we are hunting next time,
  //   otherwise we are locating next time
  hunt_or_loc_ = abs(jl-last_pos_) > kBias_ ? 0 : 1;
  last_pos_ = jl;

  // Return the lower-index edge of the stencil, without making it less than 0,
  //   and without making the upper-edge of the stencil go over [n-1]
  return std::max(0, std::min(nx-kStencilSize_, jl-((kStencilSize_-2)>>1)));
}
