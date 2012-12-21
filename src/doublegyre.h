#ifndef DOUBLEGYRE_H_
#define DOUBLEGYRE_H_

/**
 * \file doublegyre.h
 * \brief DoubleGyre velocity field
 * \author D. M. Luchtenburg
 *
 * Describes the interface of the gPC basis functions.
 */
#include "model.h"

class DoubleGyre : public Model {
public:
  DoubleGyre(double A, double eps, double omega);
  ~DoubleGyre();
  int velocity_field(double t, const double *x, double *vf) const;
  int dimen() const { return kDimen; }
private:
  const double A_;                  
  const double eps_;                   
  const double omega_;        
  static const double kPi; // DOUBLE needs to assigned in cc file!
  static const int kDimen = 2;
};

#endif  // DOUBLEGYRE_H_
