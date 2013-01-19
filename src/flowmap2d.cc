#include "flowmap2d.h"
#include "grid2d.h"
#include "array1d.h"
#include <iostream>
#include <cstdio>
#include <cassert>
#include <algorithm>
#include <cmath>


FlowMap2D::FlowMap2D(int order, double xmin, double xmax, 
                     double ymin, double ymax)
    : order_(order), xmin_(xmin), xmax_(xmax), ymin_(ymin), ymax_(ymax) {
  int num_dims = 2;
  flowmap_xdir_ = new GPCExpansion("Legendre", order_, num_dims);
  flowmap_ydir_ = new GPCExpansion("Legendre", order_, num_dims);

  // Scaling for conversion standard element
  dx_ = xmax_ - xmin_; dx_half_ = dx_ / 2.;
  dy_ = ymax_ - ymin_; dy_half_ = dy_ / 2.;
  x_mid_ = (xmax_ + xmin_) / 2.;
  y_mid_ = (ymax_ + ymin_) / 2.;
}

FlowMap2D::~FlowMap2D() {
  delete flowmap_xdir_;
  delete flowmap_ydir_;
}

void FlowMap2D::GetNodes(Grid2D &nodes) {
  int num_oned_nodes = flowmap_xdir_->num_oned_nodes();
  nodes.Resize(num_oned_nodes, num_oned_nodes);
  Array1D<double> oned_nodes(num_oned_nodes);
  oned_nodes = flowmap_xdir_->oned_nodes();
  Array1D<double> x_nodes(num_oned_nodes);
  Array1D<double> y_nodes(num_oned_nodes);
  for (int i = 0; i < num_oned_nodes; ++i) {
    StandardElToCoord(oned_nodes(i), oned_nodes(i),
                      x_nodes(i), y_nodes(i));
  }
  nodes.Mesh(x_nodes, y_nodes);
}

void FlowMap2D::SetNodalValues(const Grid2D &nodalvals) {
  int m = (int)nodalvals.Length1();
  int n = (int)nodalvals.Length2();
  // Set nodal values of flow map in (x, y)-direction
  for (int i = 0; i < m; ++i) {
    for (int j =  0; j < n; ++j) {
      int k = j + n * i;
      flowmap_xdir_->set_nodal_value(k, nodalvals.X(i,j));
      flowmap_ydir_->set_nodal_value(k, nodalvals.Y(i,j));
    }
  }
} 

void FlowMap2D::CompCoefficients() {
  // Compute coefficients
  flowmap_xdir_->CompCoefficients();
  flowmap_ydir_->CompCoefficients();
}

void FlowMap2D::Eval(Grid2D &vals) const {
  int m = vals.Length1();
  int n = vals.Length2();
  Array1D<double> point(2); // (x, y)
  double valx, valy;
  // Evaluate flow map expansion at all grid points
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      // Scale coordinates to standard element
      CoordToStandardEl(vals.X(i,j), vals.Y(i,j),
                        point(0), point(1));
      // eval flow map in (x, y)-dir
      flowmap_xdir_->EvalExpansion(point, valx);
      flowmap_ydir_->EvalExpansion(point, valy);   
      // Assign new values
      vals.X(i,j) = valx;
      vals.Y(i,j) = valy;
    }
  }
}

void FlowMap2D::CompJac(const Grid2D &vals) {
  int m = vals.Length1();
  int n = vals.Length2();
  assert((int)D11_old_.Length1() == m);
  assert((int)D11_old_.Length2() == n);
  Array1D<double> point(2); // (x, y)
  double val;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      // Scale coordinates to standard element
      CoordToStandardEl(vals.X(i,j), vals.Y(i,j),
                        point(0), point(1));
      // d/dx phi_1
      flowmap_xdir_->EvalDerExpansion(0, point, val);
      D11_new_(i,j) = val * D11_old_(i,j);
      // d/dy phi_1
      flowmap_xdir_->EvalDerExpansion(1, point, val);
      D12_new_(i,j) = val * D12_old_(i,j);
      // d/dx phi_2
      flowmap_ydir_->EvalDerExpansion(0, point, val);
      D21_new_(i,j) = val * D21_old_(i,j);
      // d/dy phi_2
      flowmap_ydir_->EvalDerExpansion(1, point, val);
      D22_new_(i,j) = val * D22_old_(i,j);
    }
  }
  // Update old vals for next iteration
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      D11_old_(i,j) = D11_new_(i,j);
      D12_old_(i,j) = D12_new_(i,j);
      D21_old_(i,j) = D21_new_(i,j);
      D22_old_(i,j) = D22_new_(i,j);
    }
  }
}


void FlowMap2D::CoordToStandardEl(const double &x, const double &y,
                                  double &xse, double &yse) const {
  xse = (x - xmin_) / dx_half_ - 1.;
  yse = (y - ymin_) / dy_half_ - 1.;
}

void FlowMap2D::StandardElToCoord(const double &xse, const double &yse,
                                  double &x, double &y) const {
  x = x_mid_ + dx_half_ * xse;
  y = y_mid_ + dy_half_ * yse;
}


void FlowMap2D::InitFTLE(const int mgrid, const int ngrid) {
  mgrid_ = mgrid;
  ngrid_ = ngrid;
  grid.Mesh(mgrid_, ngrid_, xmin_, xmax_, ymin_, ymax_);
  
  // InitJac(ngrid_, mgrid_); // ysize x xsize (mesh-grid)

  ftle_.Resize(ngrid_, mgrid_);
}


void FlowMap2D::InitJac(const int m, const int n) {
  D11_old_.Resize(m, n, 1.);
  D12_old_.Resize(m, n, 1.);
  D21_old_.Resize(m, n, 1.);
  D22_old_.Resize(m, n, 1.);
  D11_new_.Resize(m, n, 1.);
  D12_new_.Resize(m, n, 1.);
  D21_new_.Resize(m, n, 1.);
  D22_new_.Resize(m, n, 1.);
}


// void FlowMap2D::CompFTLE(const double &t,const Grid2D &vals) {
//   assert((int)vals.Length1() == ngrid_);
//   assert((int)vals.Length2() == mgrid_);
//   assert((int)ftle_.Length1() == ngrid_);
//   assert((int)ftle_.Length2() == mgrid_);
//   CompJac(vals);
//   Array2D<double> mat(2,2);
//   for (int i = 0; i < ngrid_; ++i) {
//     for (int j = 0; j < mgrid_; ++j) {
//       // Jacobian for point (i,j)
//       mat(0,0) = D11_new_(i,j); mat(0,1) = D12_new_(i,j);
//       mat(1,0) = D21_new_(i,j); mat(1,1) = D22_new_(i,j);
//       ftle_(i,j) = log(MaxSingularValue(mat)) / t;
//     }
//   }
// }


void FlowMap2D::CompFTLE(const double &t, const Grid2D &vals) {
  assert((int)vals.Length1() == ngrid_);
  assert((int)vals.Length2() == mgrid_);
  assert((int)ftle_.Length1() == ngrid_);
  assert((int)ftle_.Length2() == mgrid_);
  Array2D<double> DEL(2,2), SIG(2,2);
  double  T, D, DISCRIM, eig;
  double DURATION = t;
  for (int i = 0; i < mgrid_; ++i) {
    for (int j = 0; j < ngrid_; ++j) {
      // d/dx
      if( i== 0) {
        DEL(0,0) = (vals.X(j,i+1)-vals.X(j,i))/(grid.X(j,i+1)-grid.X(j,i));
        DEL(1,0) = (vals.Y(j,i+1)-vals.Y(j,i))/(grid.X(j,i+1)-grid.X(j,i));
      } else if (i == mgrid_ - 1) {
        DEL(0,0) = (vals.X(j,i)-vals.X(j,i-1))/(grid.X(j,i)-grid.X(j,i-1));
        DEL(1,0) = (vals.Y(j,i)-vals.Y(j,i-1))/(grid.X(j,i)-grid.X(j,i-1));
      } else {
        DEL(0,0) = (vals.X(j,i+1)-vals.X(j,i-1))/(grid.X(j,i+1)-grid.X(j,i-1));
        DEL(1,0) = (vals.Y(j,i+1)-vals.Y(j,i-1))/(grid.X(j,i+1)-grid.X(j,i-1));
      }
      // d/dy
      if (j == 0) {
        DEL(0,1) = (vals.X(j+1,i)-vals.X(j,i))/(grid.Y(j+1,i)-grid.Y(j,i));
        DEL(1,1) = (vals.Y(j+1,i)-vals.Y(j,i))/(grid.Y(j+1,i)-grid.Y(j,i));
        //
      } else if (j == ngrid_ - 1) {
        DEL(0,1) = (vals.X(j,i)-vals.X(j-1,i))/(grid.Y(j,i)-grid.Y(j-1,i));
        DEL(1,1) = (vals.Y(j,i)-vals.Y(j-1,i))/(grid.Y(j,i)-grid.Y(j-1,i));           
      } else {
        DEL(0,1) = (vals.X(j+1,i)-vals.X(j-1,i))/(grid.Y(j+1,i)-grid.Y(j-1,i));
        DEL(1,1) = (vals.Y(j+1,i)-vals.Y(j-1,i))/(grid.Y(j+1,i)-grid.Y(j-1,i));            
      }
      for (int k = 0; k < 2; ++k) {
        for (int l = 0; l < 2; ++l) {
          SIG(k,l) = DEL(0,k) * DEL(0,l) + DEL(1,k) * DEL(1,l);
        }
      }
      T = SIG(0,0) + SIG(1,1);
      D = SIG(0,0)*SIG(1,1)-SIG(0,1)*SIG(1,0);
      DISCRIM = T * T - 4 * D;
      if (DISCRIM < 0) eig = -10.;
      else eig = (T + sqrt(DISCRIM)) / 2.;
      ftle_(j,i) = (1. / DURATION) * log(sqrt(eig));
    }
  }
}




double FlowMap2D::MaxSingularValue(const Array2D<double> &mat) const {
  assert(mat.Length1() == mat.Length2());
  assert(mat.Length1() == 2);
  Array2D<double> mat2(2,2, 0.);
  Array1D<double> eigvals(2);
  // A^T A
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        mat2(i,j) += mat(k,i) * mat(k,j);
      }
    }
  }
  EigenValues(mat2, eigvals);
  return sqrt(max(eigvals(0), eigvals(1)));
}

void FlowMap2D::EigenValues(const Array2D<double> &mat, Array1D<double> &eigvals)
  const {
  assert(mat.Length1() == mat.Length2());
  assert(mat.Length1() == 2);
  assert(eigvals.Length() == 2);
  double tr, det, a, b;
  tr = mat(0,0) + mat(1,1);
  det = mat(0,0) * mat(1,1) - mat(0,1) * mat(1,0);
  a = tr; b = sqrt(tr * tr - 4. * det);
  eigvals(0) = 0.5 * (a + b);
  eigvals(1) = 0.5 * (a - b);
}


// Write Grid2D to standard output
void FlowMap2D::DumpFTLE(const char* filename) const {
  FILE* pfile = fopen(filename, "w");
  // Header
  fprintf(pfile," TITLE = \"FTLE FIELD\"\n");
  fprintf(pfile," VARIABLES = \"X\"\n");
  fprintf(pfile," \"Y\"\n");
  fprintf(pfile," \"FTLE\"\n");
  fprintf(pfile," ZONE T=\"Rectangular Zone\"\n");
  fprintf(pfile," I=%d , J=%d , ZONETYPE=Ordered\n", mgrid_, ngrid_);
  fprintf(pfile," DATAPACKING=POINT\n");
  fprintf(pfile," DT=(SINGLE SINGLE SINGLE)\n");
  // Data
  for (int i = 0; i < (int)ngrid_; i++) {
    for (int j = 0; j < (int)mgrid_; j++) {
      fprintf(pfile, "%.5e %.5e %.5e \n", grid.X(i,j), grid.Y(i,j), ftle_(i,j));
    }
  }
  fclose(pfile);
}
