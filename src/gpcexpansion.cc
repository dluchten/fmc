#include "gpcexpansion.h"
#include "legendrebasis.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include "grid2d.h"


/**
 * \todo maybe implement hash table
 *
 */
GPCExpansion::GPCExpansion(string type, int order, int num_dims,
                           int num_oned_nodes) 
    : order_(order), num_oned_nodes_(num_oned_nodes), num_dims_(num_dims) {
  // Select GPC Basis type. Hermite not supported yet.
  if (type.compare("Legendre") == 0) {
    basis_ = new LegendreBasis(order, num_oned_nodes);
  } 
  else if (type.compare("Hermite") == 0) {
    std::cout << "Hermite not supported yet." << endl;
    exit (EXIT_FAILURE);
  }
  else {
    std::cout << "Error in GPCExpansion: Invalid choice for type." << endl;
    exit(EXIT_FAILURE);
  }
  num_oned_basisfncts_ = order_ + 1;
  num_nodes_ = (int) pow((double)num_oned_nodes_, num_dims_);

  // std::cout << "CHECK - num_nodes_ = " << num_nodes_ << endl;
  
  nodal_values_.Resize(num_nodes_, 0.);

  CompMultiIndicesExp();
  CompMultiIndicesNodes();
  InitProjection();

  // Array for storage 
  // GPC Basis types. Hermite not supported yet.
  // enum {Legendre, Hermite};
  // switch (type) {
  //   case Legendre:
  //     basis_ = new LegendreBasis(order, num_nodes);
  //     break;

  //   case Hermite:
  //     std::cout << "Hermite not supported yet." << endl;
  //     exit (EXIT_FAILURE);
  //   default:
  //     std::cout << "Error in GPCExpansion: Basis type not supported." << endl;
  //     exit (EXIT_FAILURE);
  // }
}

GPCExpansion::GPCExpansion(string type, int order, int num_dims)
  : order_(order),  num_dims_(num_dims) {

  // Select GPC Basis type. Hermite not supported yet.
  if (type.compare("Legendre") == 0) {
    basis_ = new LegendreBasis(order);
  } 
  else if (type.compare("Hermite") == 0) {
    std::cout << "Hermite not supported yet." << endl;
    exit (EXIT_FAILURE);
  }
  else {
    std::cout << "Error in GPCExpansion: Invalid choice for type." << endl;
    exit(EXIT_FAILURE);
  }
  num_oned_basisfncts_ = order_ + 1;
  num_oned_nodes_ = basis_->num_nodes();
  num_nodes_ = (int) pow((double)num_oned_nodes_, num_dims_);

  // std::cout << "CHECK - num_nodes_ = " << num_nodes_ << endl;
  

  nodal_values_.Resize(num_nodes_, 0.);

  CompMultiIndicesExp();
  CompMultiIndicesNodes();
  InitProjection();
}


GPCExpansion::~GPCExpansion() {
  multi_index_exp_.Clear(); 
  multi_index_nodes_.Clear(); 
  omega_.Clear(); 
  nodal_values_.Clear(); 
  delete basis_;
  coefficients_.Clear(); 
}
                               
void GPCExpansion::ComputeDimExp() {
  int num = 1, den = 1;
  int m = min(order_, num_dims_);
  for( int i = 0; i < m; ++i) {
      num = num * (order_ + num_dims_ - i);
      den = den * (i + 1);
  }
  dim_exp_ = num / den;
}

void GPCExpansion::CompMultiIndicesExp() {
  // Compute number of expansion coefficients: (order + num_dim) over (order) 
  ComputeDimExp();
  // Init multi-index with zeros
  multi_index_exp_.Resize(dim_exp_, num_dims_, 0);
  int index = 0, isum = 0;
  // Work arrays
  Array1D<int> ic(num_dims_,0);
  Array1D<int> ict(num_dims_,0);
  
  // Compute multi-index
  // 0th-order: (0,...,0) 
  index = 0;
  for (int idim = 0; idim < num_dims_; ++idim) {
    multi_index_exp_(index,idim) = 0;
  }
  if (order_ > 0 ) {
  // 1st-order: identity matrix
    for (int idim = 0; idim < num_dims_; idim++) {
      index++;
      multi_index_exp_(index, idim) = 1;
      ic(idim) = 1;
    }
  }
  // higher order terms
  if (order_ > 1) {
    for (int iord = 2; iord < order_+1; iord++) {
      int lessiord = index; // number of terms of order less than iord
      for (int idim = 0; idim < num_dims_; idim++) {
        isum = 0;
        for(int ii = idim; ii < num_dims_; ii++) {
          isum += ic(ii);
        }
        ict(idim) = isum;
      }
      // Update 
      for (int idim = 0; idim < num_dims_; idim++) {
        ic(idim)=ict(idim);
      }

      // std::cout << "ic--------" << "\n";
      // for (int j = 0; j < num_dims_; j++) {
      //   std::cout << ic(j);
      // }
      // std::cout << "\n";

      for (int idimm = 0; idimm <num_dims_; idimm++) {
        for (int ii = lessiord -ic(idimm) + 1; ii < lessiord+1; ii++){
          index++;
          for (int idim = 0; idim < num_dims_; idim++) {
            multi_index_exp_(index,idim) = multi_index_exp_(ii,idim);

            // std::cout << "A multi-index(" << index+1 << "," << idim+1 << ") = " << multi_index_(index,idim) << "\n";

          }
          multi_index_exp_(index,idimm) = multi_index_exp_(index,idimm) + 1;

            // std::cout << "B multi-index(" << index+1 << "," << idimm+1 << ") = " << multi_index_(index,idimm) << "\n";

        }
      }
    }
  }
  // Check if index equals dimension of expansion
  assert(index + 1 == dim_exp_);

  // std::cout << "Multi-indices" << "\n";
  // std::cout << "dimension of exp. = " << dim_exp_ << "\n";
  // for (int i = 0; i < dim_exp_; ++i) {
  //   std::cout << i << "   ";
  //   for (int j = 0; j < num_dims_; ++j) {
  //     std::cout << multi_index_exp_(i,j) << " ";
  //   }
  //   std::cout << "\n";
  // }

}

void GPCExpansion::CompMultiIndicesNodes() {
  // multi-index quad. nodes
  multi_index_nodes_.Resize(num_nodes_, num_dims_, 0); 
  // work array
  Array1D<int> ic(num_nodes_, 0); 

  int M, N; M = N = 1;
  for (int i = 0; i < num_dims_; ++i) {
    if ( i > 0 ) {
      N *= num_oned_nodes_; 
      M = N / num_oned_nodes_;
    }
    for (int j = 0; j < num_nodes_; ++j) {
      if (i == 0) {
        ic(j) = j;
      } else {
        ic(j) = ic(j) - multi_index_nodes_(j,i-1) * M;
      }
      multi_index_nodes_(j,i) = (ic(j) / N) % num_oned_nodes_;
    }
  }

  // std::cout << "num_oned_nodes_" << num_oned_nodes_ <<"\n";
  // std::cout << "indices" << "\n";
  // int k = 0;
  // for (int i = 0; i < num_nodes_; ++i) {
  //   std::cout << "k = " << k << ", ";
  //   k++;
  //   for (int j = 0; j < num_dims_; ++j) {
  //     std::cout << multi_index_nodes_(i,j);
  //   }
  //   std::cout << "\n";
  // }

}


void GPCExpansion::set_nodal_value(int index, double value) {
  assert(index <= num_nodes_);
  nodal_values_(index) = value;
}

void GPCExpansion::InitProjection() {
  omega_.Resize(dim_exp_, num_nodes_, 1.);
  int index_node, index_exp;
  for (int i = 0; i < dim_exp_; ++i) {
    for (int k = 0; k < num_nodes_; ++k) {
      // Loop over all dimensions for expansion and node index
      for (int idim = 0; idim < num_dims_; ++idim) {
        index_exp = multi_index_exp_(i, idim);
        index_node = multi_index_nodes_(k, idim);
        
        // Checks can be commented out
        assert(index_exp <= order_);
        assert(index_node <= num_oned_nodes_);
        // Checks can be commented out

        omega_(i,k) *= basis_->omega(index_node, index_exp);
      }
    }
  }
}

void GPCExpansion::CompCoefficients() {
  coefficients_.Resize(dim_exp_, 0.);
  for (int i = 0; i < dim_exp_; ++i) {
    // summation, discrete Galerkin projection
    for (int k = 0; k < num_nodes_; ++k) {
      coefficients_(i) += omega_(i,k) * nodal_values_(k);
    }
  }

  // std:: cout << "nodal values" << "\n";
  // for (int k = 0; k < num_nodes_; ++k) {
  //   std::cout << "nodal_values_(" << k << ") = " << nodal_values_(k) << "\n";
  // }
  // std::cout << "\n";
}

void GPCExpansion::PrintCoefficients() {
  std:: cout << "GPCExpansion coefficients" << "\n";
  for (int i = 0; i < dim_exp_; ++i) {
    printf("% 6.4f \n", coefficients_(i));
  }
}

void GPCExpansion::EvalExpansion(const Array1D<double> &point, double &value) {
  // Check if length point equals number of dimensions
  assert(point.Length() == (unsigned)num_dims_);
  // Evaluate 1D basis functions at point values
  Array2D<double> basis_evals_alldim(num_oned_basisfncts_, num_dims_, 0.);
  Array1D<double> basis_evals(num_oned_basisfncts_, 0.);
  //  phi_i(z1), phi_i(z2), 
  for (int i = 0; i < num_dims_; ++i) {
    basis_->EvalBasisFunctions(point(i), basis_evals);
    // Assign values to array
    for (int j = 0; j < num_oned_basisfncts_; ++j) {
      basis_evals_alldim(j, i) =  basis_evals(j);
    }
  }
  // Sum all terms: sum_I coefficients_I * PHI_I(point)
  value = 0.;
  for (int i = 0; i < dim_exp_; ++i) {
    // get basis function PHI_I(Z) = phi_i1(z1) * phi_i2(z2) * .. 
    double PHI = 1.;
    for(int idim = 0; idim < num_dims_; ++idim) {
      // index basis function of dimension idim
      int basis_index = multi_index_exp_(i, idim);
      PHI *= basis_evals_alldim(basis_index ,idim);
    }
    value += coefficients_(i) * PHI;
  }
}


void GPCExpansion::EvalDerExpansion(const int idim, const Array1D<double> &point, double &value) {
  assert(point.Length() == (unsigned)num_dims_);
  assert(idim < num_dims_);

  // Evaluate 1D basis functions at point values
  Array2D<double> basis_evals_alldim(num_oned_basisfncts_, num_dims_, 0.);
  Array1D<double> basis_evals(num_oned_basisfncts_, 0.);
  //  phi_i(z1), phi_i(z2), d / zi (phi(zi)), ...
  for (int i = 0; i < num_dims_; ++i) {
    if( i == idim ) {
      basis_->EvalDerBasisFunctions(point(i), basis_evals);
    } else {
      basis_->EvalBasisFunctions(point(i), basis_evals);
    }
    // Assign values to array
    for (int j = 0; j < num_oned_basisfncts_; ++j) {
      basis_evals_alldim(j, i) =  basis_evals(j);
    }
  }
  // Sum all terms: sum_I coefficients_I * PHI_I(point)
  value = 0.;
  for (int i = 0; i < dim_exp_; ++i) {
    // get basis function PHI_I(Z) = phi_i1(z1) * phi_i2(z2) * .. 
    double PHI = 1.;
    for(int idim = 0; idim < num_dims_; ++idim) {
      // index basis function of dimension idim
      int basis_index = multi_index_exp_(i, idim);
      PHI *= basis_evals_alldim(basis_index ,idim);
    }
    value += coefficients_(i) * PHI;
  }
}
