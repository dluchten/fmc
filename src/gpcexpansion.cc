#include "gpcexpansion.h"
#include "legendrebasis.h"
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>



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
  num_nodes_ = pow(num_oned_nodes_, num_dims_);
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
  num_oned_nodes_ = basis_->num_nodes();
  num_nodes_ = pow(num_oned_nodes_, num_dims_);
  nodal_values_.Resize(num_nodes_, 0.);

  CompMultiIndicesExp();
  CompMultiIndicesNodes();
  InitProjection();

  ///////////////////////
 CompCoeffients();
  ///////////////////////
}


GPCExpansion::~GPCExpansion()
{}
                               


void GPCExpansion::CompMultiIndicesExp() {
  // Compute number of expansion coefficients: (order + num_dim) over (order) 
  int num = 1, den = 1;
  for (int i = 1; i < order_ + 1; ++i) {
    num *= num_dims_ + i;
    den *= i;
  } 
  dim_exp_ = num / den;
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

 
  // std::cout << "dimension of exp. = " << dim_exp_ << "\n";
  // for (int i = 0; i < dim_exp_; ++i) {
  //   std::cout << i << " ";
  //       for (int j = 0; j < num_dims_; ++j) {
  //         std::cout << multi_index_(i,j);
  //       }
  //       std::cout << "\n";
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
  //   std::cout << k << " ";
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
      // Loop over all dimensions for expansion index
      for (int idim = 0; idim < num_dims_; ++idim) {
        index_exp = multi_index_exp_(i, idim);
        // Loop over all nodes for node index
        for (int inode = 0; inode < num_oned_nodes_; ++inode) {
          index_node = multi_index_nodes_(k, inode);
          omega_(i,k) *= omega_(i,k) *= basis_->omega(index_node, index_exp);
        }
      }
    }
  }
}

void GPCExpansion::CompCoeffients() {
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

  // std:: cout << "coefficients" << "\n";
  // for (int i = 0; i < dim_exp_; ++i) {
  //   std::cout << "coefficients_(" << i << ") = " << coefficients_(i) << "\n";
  // }



}
