#ifndef GRID2D_H_
#define GRID2D_H_

/**
 * \file grid2d.h
 * \brief Class for storage of 2D grid (matrix) 
 * \author D. M. Luchtenburg
 *
 * Implementation of 2D grid to store data.
 */
#include <vector>
#include <cassert>
#include "array1d.h"
#include "array2d.h"

using namespace std;

/**
 * \brief Stores data in 2D array (matrix) of specified type T
 */
class Grid2D {
 public:
  /// \brief Default constructor
  Grid2D() : m_(0), n_(0) {}

  /// \brief Constructor that allocates the memory
  Grid2D(const size_t &m, const size_t &n) : m_(m), n_(n) {
    x_.Resize(m_, n_);
    y_.Resize(m_, n_);
  }

  /// \brief Constructor that sets uniformly sampled mesh-grid
Grid2D(const size_t &m, const size_t &n, 
       const double &xmin, const double &xmax, 
       const double &ymin, const double &ymax);

  /** \brief Constructor that allocates memory, and initializes data to mesh-grid values specified by x and y
   */
  Grid2D(const Array1D<double> &x, const Array1D<double> &y);

  /// \brief Destructor
  ~Grid2D() { x_.Clear(); y_.Clear(); }

  /// \brief Mesh-grid based on values specified by x and y
  void Mesh(const size_t &m, const size_t &n, 
            const double &xmin, const double &xmax, 
            const double &ymin, const double &ymax);
  void Mesh(const Array1D<double> &x, const Array1D<double> &y);

  /* 
   * \brief Resizes Grid2D
   */
  void Resize(const size_t &m, const size_t &n) {
    m_ = m;
    n_ = n;
    x_.Resize(m_, n_);
    y_.Resize(m_, n_);    
  }

  /// \brief Returns size in the 1st-direction
  size_t Length1() const { return m_; }

  /// \brief Returns size in the 2nd-direction
  size_t Length2() const { return n_; }

  /// \brief Returns size
  size_t Length() const { return m_ * n_; }
 
  /**
   * \brief Assign / access Array2D values (row-major order)
   */
  double& X(size_t i, size_t j) { return x_(i,j); }
  const double& X(size_t i,size_t j) const { return x_(i,j); }
  Array2D<double>& X() { return x_; }
  const Array2D<double>& X() const { return x_; }
  double& Y(size_t i, size_t j) { return y_(i,j); }
  const double& Y(size_t i,size_t j) const { return y_(i,j); }
  Array2D<double>& Y() { return y_; }
  const Array2D<double>& Y() const { return y_; }

  Grid2D& operator=(const Grid2D &rhs);

  // Write grid to ASCII file
  void Dump(const char* filename) const;

  private:
    size_t m_; ///> Number of elements in 1st-dimension
    size_t n_; ///> Number of elements in 2nd-dimension

    /**
     * \brief Data in array with size = m_ * n_
     *
     * The data is stored in row-major order (2nd index is fast-varying).
     */
    Array2D<double> x_;
    Array2D<double> y_;
};

#endif // GRID2D_H_
