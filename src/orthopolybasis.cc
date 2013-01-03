#include "orthopolybasis.h"
#include <string>
#include "array2d.h"
#include <iostream>


OrthoPolyBasis::OrthoPolyBasis(string type, int order)
    : type_(type), order_(order) {
  num_basis_fncts_ = order_ + 1;
}

OrthoPolyBasis::OrthoPolyBasis(string type, int order, int num_nodes)
    : type_(type), order_(order), num_nodes_(num_nodes) {
  num_basis_fncts_ = order_ + 1;
}

void OrthoPolyBasis::CompOmegaGP() {
  omega_.Resize(num_nodes_, num_basis_fncts_, 0.);
  for (int i = 0; i < num_nodes_; ++i) {
    for (int j = 0; j < num_basis_fncts_; ++j) {
      omega_(i,j) = psi_(i,j) * weights_(i) / gamma_(j);
    }
  }

  // std:: cout << "omega_oned" << "\n";
  // for (int i = 0; i < num_nodes_; ++i) {
  //   for (int j = 0; j < order_; ++j) {
  //     std::cout << omega_(i,j) << " ";
  //   }
  //     std:: cout << "\n";
  // }
  //     std:: cout << "\n";

}



