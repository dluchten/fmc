#ifndef _GPCBASIS_H_
#define _GPCBASIS_H_

/**
 * \file gpcbasis.h
 * \brief GPC basis class (Generalized Polynomial Chaos)
 * \author D. M. Luchtenburg
 *
 * Describes the interface of the gPC basis functions.
 */
#include <string>
#include "array1d.h"
#include "array2d.h"

using namespace std;

class GPCBasis {
public:
  /// \brief Constructor
  GPCBasis(string type, int order);
  /// \brief Destructor
  virtual ~GPCBasis() {}
  /// \brief Evaluate all basis functions for a given value z.
  virtual void EvalBasisFunctions(const double &z,
                                  Array1D<double> &basis_evals) const = 0;
  /// \brief Get number of quadrature nodes.
  int num_nodes() const { return num_nodes_; }
  int order() const { return order_; }
  /// \brief Get quadrature node #i.
  const double& nodes(int i) { return nodes_(i); }
  /// \brief Get quadrature weight #i.
  const double& weights(int i) { return weights_(i); }
protected: // give derived classes access
  /// \brief Initialize quadrature nodes and weights.
  virtual void InitNodesWeights() = 0;
  /// \brief Evaluate basis function at quadrature nodes.
  virtual void EvalBasisFunctionsAtNodes() = 0;
  /// \brief Compute squared norms of basis functions.
  virtual void CompBasisFunctionNormsSquared() = 0;
  const string type_; ///< type of basis functions
  const int order_; ///< order of basis
  int num_nodes_; ///< number of quadrature nodes
  Array1D<double> nodes_; ///< quadrature nodes
  Array1D<double> weights_; ///< quadrature weights
  Array1D<double> gamma_; ///< basis function norms squared
  Array2D<double> psi_; ///< psi_(i,j) basis function of degree j at node i
};


#endif  //  _GPCBASIS_H_ 
