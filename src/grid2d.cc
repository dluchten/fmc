#include "grid2d.h"
#include <cstdio>

Grid2D::Grid2D(const Array1D<double> &x, const Array1D<double> &y) {
  m_ = y.Length();
  n_ = x.Length();
  x_.Resize(m_, n_);
  y_.Resize(m_, n_);
  for (size_t i = 0; i < m_; ++i) {
    for (size_t j = 0; j < n_; ++j) {
      x_(i,j) = x(j);
      y_(i,j) = y(i);
    }
  }
}


Grid2D::Grid2D(const size_t &m, const size_t &n, 
               const double &xmin, const double &xmax, 
               const double &ymin, const double &ymax) : m_(n), n_(m) {
  x_.Resize(m_, n_);
  y_.Resize(m_, n_);
  double dx = (xmax - xmin) / (n_ - 1);
  double dy = (ymax - ymin) / (m_ - 1);
  for (size_t i = 0; i < m_; ++i) {
    for (size_t j = 0; j < n_; ++j) {
      x_(i,j) = xmin + j * dx;
      y_(i,j) = ymin + i * dy;
    }
  }
}

void Grid2D::Mesh(const Array1D<double> &x, const Array1D<double> &y) {
  assert(m_ == y.Length());
  assert(n_ == x.Length());
  for (size_t i = 0; i < m_; ++i) {
    for (size_t j = 0; j < n_; ++j) {
      x_(i,j) = x(j);
      y_(i,j) = y(i);
    }
  }
}

void Grid2D::Mesh(const size_t &m, const size_t &n, 
                  const double &xmin, const double &xmax, 
                  const double &ymin, const double &ymax) {
  m_ = n;
  n_ = m;
  x_.Resize(m_, n_);
  y_.Resize(m_, n_);
  double dx = (xmax - xmin) / (n_ - 1);
  double dy = (ymax - ymin) / (m_ - 1);
  for (size_t i = 0; i < m_; ++i) {
    for (size_t j = 0; j < n_; ++j) {
      x_(i,j) = xmin + j * dx;
      y_(i,j) = ymin + i * dy;
    }
  }
}



// Write Grid2D to standard output
void Grid2D::Dump(const char* filename) const {
  FILE* pfile = fopen(filename, "w");
  for (int i = 0; i < (int)m_; i++) {
    for (int j = 0; j < (int)n_; j++) {
      fprintf(pfile, "%.5e %.5e  \n", x_(i,j), y_(i,j));
    }
  }
  fclose(pfile);
}
