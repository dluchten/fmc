#ifndef _LEGENDREBASIS_H_
#define _LEGENDREBASIS_H_

/** 
 * \file legendrebasis.h
 * \brief Legendre basis function class (for uniform distribution)
 * \author D. M. Luchtenburg
 *
 * This class describes the implementation of the Legendre basis functions.
 * \note We focus on non-intrusive GPC (stochastic collocation).  
 * Number of quadrature nodes equals: order expansion + 1.
 * \todo Make this more flexible in future 
 * (e.g. higher order integrands)
 */
#include "gpcbasis.h"
#include "array1d.h"

class LegendreBasis : public GPCBasis {
 public:
  LegendreBasis(int order);
  ~LegendreBasis();
  /// \brief Evaluates all basis functions at z.
  void EvalBasisFunctions(const double &z, Array1D<double> &basis_evals) const;
 private:
  /*
   * \brief Initialize quadrature nodes and weights 
   * for Gauss-Lobatto-Legendre integration.
   */
  void InitNodesWeights();
  /// \brief Evaluate basis function at quadrature nodes.
  void EvalBasisFunctionsAtNodes();
  /// \brief Compute squared norms of basis functions.
  void CompBasisFunctionNormsSquared(); 
};

#endif  // _LEGENDREBASIS_H_

