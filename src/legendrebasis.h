#ifndef _LEGENDREBASIS_H_
#define _LEGENDREBASIS_H_

/** 
 * \file legendrebasis.h
 * \brief Legendre basis function class
 * \author D. M. Luchtenburg
 *
 * This class describes the implementation of the Legendre basis functions.
 * Includes: basis functions and derivatives, quadrature rule (Gauss-Lobatto)
 */
#include "orthopolybasis.h"
#include "array1d.h"

class LegendreBasis : public OrthoPolyBasis {
 public:
  LegendreBasis(int order);
  LegendreBasis(int order, int num_nodes);
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

