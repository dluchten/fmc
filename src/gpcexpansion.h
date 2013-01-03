#ifndef _GPCEXPANSION_H_
#define _GPCEXPANSION_H_

#include <string>
#include "array1d.h"
#include "array2d.h"
#include "orthopolybasis.h"

class GPCExpansion {
public:
  /// \brief Constructors
  GPCExpansion(string type, int order, int num_dim, int num_oned_nodes);
  GPCExpansion(string type, int order, int num_dim);
  /// \brief Destructor
  ~GPCExpansion();
  /// \brief Evaluate all basis functions for a given value z.PCExpansion(string type, int order, int num_dim, int num_oned_nodes);

  /// \brief Get number of quadrature nodes.
  int num_nodes() const { return basis_->num_nodes(); }
  /// \brief Get order (maximum degree of basis function)
  int order() const { return  basis_->order(); }

  /** \brief Set value at quadrature node with multi-index (i1, i2, i3,...)
   *
   * index = i1 + i2 * N1 + i3 * N1 * N2 + ...
   */
  void set_nodal_value(int index, double value);


  // if all are computed: compute expansion coefficients

  /// \brief Compute expansion coefficients
  void CompCoeffients();


  /// \brief EvalExpansion()
  /// \brief get multi index

  /// \brief get number of expansion coefficients

  // Get number of quadrature nodes.



private:  
  /// \brief Compute multi-indices for GPC expansion
  void CompMultiIndicesExp(); 
  /// \brief Compute multi-indices for quadrature nodes
  void CompMultiIndicesNodes();
  /// \brief Initialize weights for discrete Galerkin projection
  void InitProjection();

  const int order_; ///< order
  int num_oned_nodes_; ///< number of quadrature nodes in one dimension
  int num_nodes_; ///< total number of quadrature nodes
  const int num_dims_; ///< number of dimensions
  int dim_exp_; ///<  dimension of expansion (number of coefficients)

  Array2D<int> multi_index_exp_; ///< multi-indices for GPC exp.
  Array2D<int> multi_index_nodes_; ///< multi-indices for quad. nodes
  /// \brief weights for discrete Galerkin projection
  Array2D<double> omega_;   // omega(exp_index,node_index)
  Array1D<double> nodal_values_; ///< function values at nodes
  OrthoPolyBasis *basis_;   ///< Basis functions (same in every dimension)
  Array1D<double> coefficients_; ///< expansion coefficients
};

#endif  //  _GPCEXPANSION_H_
