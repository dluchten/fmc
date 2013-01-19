#ifndef DOUBLEGYRE_H_
#define DOUBLEGYRE_H_

/**
 * \file doublegyre.h
 * \brief DoubleGyre velocity field
 * \author D. M. Luchtenburg
 *
 */
#include "velocityfield.h"
#include <cmath>

const double kPi = acos(-1.);

class DoubleGyre : public VelocityField {
public:
  /// \brief Default constructor
  DoubleGyre(double A = 0.1, double eps = 0.1, double omega = 2 * kPi / 10.);
  /// \brief destructor
  ~DoubleGyre();
  /// \brief Get velocity field vf(t,x)
  void Get(double t, const double *x, double *vf) const;
  int dimen() const { return kDimen; }
private:
  const double A_;                  
  const double eps_;                   
  const double omega_;        
  static const int kDimen = 2;
};

#endif  // DOUBLEGYRE_H_
