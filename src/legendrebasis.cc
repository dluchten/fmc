#include "array1d.h"
#include "array2d.h"
#include "legendrebasis.h"
#include <string>
#include <cmath>
#include <cassert>

#include <cstdio>
#include <iostream>


/// \todo FIX ANNOYING: order_ + 1 ---> num_basis_functions
LegendreBasis::LegendreBasis(int order, int num_nodes)
    : OrthoPolyBasis("Legendre", order, num_nodes)
{
  // assert(num_nodes >= order_ + 2);
  InitNodesWeights();
  EvalBasisFunctionsAtNodes();
  CompBasisFunctionNormsSquared();
  CompOmegaGP();
}

LegendreBasis::LegendreBasis(int order)
  : OrthoPolyBasis("Legendre", order)
{
  // Set number of nodes (interpolating polynomial through the GLL nodes).
  num_nodes_ = order_ + 1;
  InitNodesWeights();
  EvalBasisFunctionsAtNodes();
  CompBasisFunctionNormsSquared();
  CompOmegaGP();
}

LegendreBasis::~LegendreBasis()
{}

/*
 * \brief Evaluates all basis functions at z.
 *
 * \note Might be better to call this class LegendreLobatto since end points
 * are included.
 *
 * Use recursion
 * P_n(x) = (2*n - 1) / n * x * P_{n-1}(x) - (n-1) / n * P_{n-2}(x), n > 1,
 * with start values: P_0 = 1, P_1 = x
 */
void LegendreBasis::EvalBasisFunctions(const double &z,
                                       Array1D<double> &basis_evals) const {
  // Check size
  assert(basis_evals.Length() == (unsigned)(num_basis_fncts_));
  // Use recursion
  basis_evals(0) = 1.;
  if (order_ > 0) {
    basis_evals(1) = z;
    for (int i = 2; i < order_ + 1; ++i) {
      basis_evals(i) = ((1. - i) * basis_evals(i - 2) +
                        (2. * i - 1) * z * basis_evals(i - 1)) / i;
    }
  }
}

/*
 * \brief Evaluates all derivatives of basis functions at z.
 *
 * Use recursion (see also EvalBasisFunctions above)
 * (1 - x^2) * P'_n(x) = n P_{n-1}(x) - n * x * P_{n}(x), n > 1,
 */
void LegendreBasis::EvalDerBasisFunctions(const double &z,
                                          Array1D<double> &basis_evals) const {
  // Evaluate all Legendre polynomials at z
  Array1D<double> basis_evals_orig(num_basis_fncts_);
  EvalBasisFunctions(z, basis_evals_orig);
  // Check size
  assert(basis_evals.Length() == (unsigned)(num_basis_fncts_));  
  // Compute derivatives
  basis_evals(0) = 0.;
  if (order_ > 0) {
    if (z == -1.) {
      int sign = 1;
      for (int i = 1; i < order_ + 1; ++i) {
        basis_evals(i) = sign * i * (i +  1) / 2.;
        sign *= -1;
      }
    } else if (z == 1.) {
      for (int i = 1; i < order_ + 1; ++i) {
        basis_evals(i) = i * (i +  1) / 2.;
      }
    } else {
      for (int i = 1; i < order_ + 1; ++i) {
        basis_evals(i) = (1. / (1. - z * z)) *
          (i * basis_evals_orig(i-1) - i * z * basis_evals_orig(i));
      }
    }
  }
}


/*
 * \brief Compute Gauss-Lobatto-Legendre quadrature (nodes, weights)
 *
 * The Gauss-Lobatto-Legendre quadrature nodes are the zeros of
 * (1 - x^2) * P'_{n - 1}(x), where P_n is a n-th degree Legendre
 * polynomial.
 * \todo Write function to get machine precision, and use it to set tolerance
 */
void LegendreBasis::InitNodesWeights() {
  static const double kPi = abs(acos(-1.)); // pi
  nodes_.Resize(num_nodes_, 0.);
  weights_.Resize(num_nodes_);
  const double kTol = 1.e-15; // tolerance
  int m = (num_nodes_ + 1) / 2; // exploit symmetry of roots
  double z, zOld; // store new and old guesses zeros
  double p1, p2, p3; // store evalations of Legendre polynomials
  
  // Find the roots of P'_{n - 1}(x) using Newton's method,
  for (int i = 0; i < m; ++i) {
    // Initial guess, use Gauss-Lobatto-Chebyshev nodes
    z = cos(i * kPi / (num_nodes_ - 1));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (int j = 0; j < num_nodes_ - 1; ++j) {
        p3 = p2;
        p2 = p1;
        p1 = ( - j * p3 + (2. * j + 1) * z * p2 ) / (j + 1); // P_j(z)
      }
      // p1'  =  (n - 1) * (z * p1 - p2) / (z^2 - 1) 
      // p1'' =  ((n - 1) * n * p1 - 2 * z * p1') / (z^2 - 1)  (Legendre DE)
      // neglect 2nd term in numerator since p1' -> 0 and abs(z) <= 1
      // Newton step
      zOld = z;
      z = zOld - (z * p1 - p2) / (num_nodes_ * p1);
    } while (abs(z - zOld) > kTol);
    // Quadrature nodes
    nodes_(i) = -z;
    nodes_(num_nodes_ - i - 1) = z; // use symmetry
    // Quadrature weights
    weights_(i) = 2. / ((num_nodes_ - 1) * num_nodes_ * p1 * p1);
    weights_(num_nodes_ - i - 1) = weights_(i); // symmetry of roots
  }
  // Set Lobatto (boundary) nodes exactly
  nodes_(0) = -1.;
  nodes_(num_nodes_ - 1) = 1.;
}

void LegendreBasis::EvalBasisFunctionsAtNodes() {
  psi_.Resize(num_nodes_, order_ + 1, 0.);
  Array1D<double> basis_evals(order_ + 1);
  // Evaluate all basis functions for each node i
  for (int i = 0; i < num_nodes_; ++i) {
    EvalBasisFunctions(nodes_(i), basis_evals);
    // Add all basis function evaluations of node i to psi_
    for (int j = 0; j < order_ + 1; j++) {
      psi_(i,j) = basis_evals(j);
    }
  } 
}

/*
 * \brief Compute squared norms of basis functions, 
 * using quadrature.
 */
void LegendreBasis::CompBasisFunctionNormsSquared() {
  gamma_.Resize(order_ + 1);
  double tmp;
  // Compute squared norms for all basis functions i
  for (int i = 0; i < order_ + 1; ++i) {
    // Sum over all nodes
    tmp = 0.;
    for (int j = 0; j < num_nodes_; ++j) {
      tmp += psi_(j, i) * psi_(j, i) * weights_(j);
    }
    gamma_(i) = tmp;
  }
}
