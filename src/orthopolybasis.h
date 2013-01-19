#ifndef _ORTHOPOLYBASIS_H_
#define _ORTHOPOLYBASIS_H_

/**
 * \file orthobasis.h
 * \brief Orthogonal Polynomial Basis
 * \author D. M. Luchtenburg
 *
 * Describes the interface of orthogonal polynomial basis functions.
 */
#include <string>
#include "array1d.h"
#include "array2d.h"
#include <iostream>
using namespace std;

class OrthoPolyBasis {
public:
  /// \brief Constructors
  OrthoPolyBasis(string type, int order);
  OrthoPolyBasis(string type, int order, int num_nodes);
  /// \brief Destructor
  virtual ~OrthoPolyBasis() {}
  /// \brief Evaluate all basis functions for a given value z.
  virtual void EvalBasisFunctions(const double &z, Array1D<double> &basis_evals)
      const = 0;
  /// \brief Evaluates all derivatives of basis functions at z.
  virtual void EvalDerBasisFunctions(const double &z, 
                                     Array1D<double> &basis_evals) const = 0;
  /// \brief Get number of quadrature nodes.
  int num_nodes() const { return num_nodes_; }
  /// \brief Get order (maximum degree of basis function)
  int order() const { return order_; }
  /// \brief Get quadrature node #i.
  const double& nodes(int i) const { return nodes_(i); }
  /// \brief Get all quadrature nodes
  const Array1D<double>& nodes() const { return nodes_; }

  /// \brief Get quadrature weight #i.
  const double& weights(int i) const { return weights_(i); }
  /// \brief Get basis function j at node node i                                
  const double& basisfunctionsatnodes(const int i, const int j) const
                                  { return psi_(i,j); }
  const double& omega(int i, int j) const { return omega_(i,j); }
protected: // give derived classes access
  /// \brief Initialize quadrature nodes and weights.
  virtual void InitNodesWeights() = 0;
  /// \brief Evaluate basis function at quadrature nodes.
  virtual void EvalBasisFunctionsAtNodes() = 0;
  /// \brief Compute squared norms of basis functions.
  virtual void CompBasisFunctionNormsSquared() = 0;
  /// \brief Compute omega for discrete Galerkin projection
  void CompOmegaGP();

  const string type_; ///< type of basis functi
  const int order_; ///< order of basisons
  int num_basis_fncts_; ///< number of basis functions

  int num_nodes_; ///< number of quadrature nodes
  Array1D<double> nodes_; ///< quadrature nodes
  Array1D<double> weights_; ///< quadrature weights
  Array1D<double> gamma_; ///< basis function norms squared
  Array2D<double> psi_; ///< psi_(i,j) basis function of degree j at node i

  /// \brief omega_(node i, basis j) = phi_j(i) * w(i) / gamma(j)
  Array2D<double> omega_; 
  
};


#endif  //  _ORTHOPOLYBASIS_H_ 
