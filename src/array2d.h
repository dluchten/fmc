#ifndef ARRAY2D_H_
#define ARRAY2D_H_

/**
 * \file array2d.h
 * \brief Class for storage of 2D array (matrix) of specified type T
 * \author D. M. Luchtenburg
 *
 * Implementation of 2D array to store data.
 */
#include <vector>
#include <cstddef> 
using namespace std;

/**
 * \brief Stores data in 2D array (matrix) of specified type T
 */
template <typename T>
class Array2D {
 public:
  /// \brief Default constructor
  Array2D() : m_(0), n_(0) {}

  /// \brief Constructor that allocates the memory
  Array2D(const size_t &m, const size_t &n) : m_(m), n_(n) {
    data_.resize(m_ * n_);
  }

  /** \brief Constructor that allocates memory, and initializes data to constant
   * value
   */
  Array2D(const size_t &m, const size_t &n, const T &t) : m_(m), n_(n) {
    data_.resize(m_ * n_, t);
  }

  /// \brief Destructor
  ~Array2D() { data_.clear(); }

  /* /// \brief Function to clear the memory */
  /* void Clear() { m_ = 0; n_ = 0; data_.clear(); } */

  /// \brief Returns size in the 1st-direction
  size_t Length1() const { return m_; }

  /// \brief Returns size in the 2nd-direction
  size_t Length2() const { return n_; }

  /// \brief Returns size
  size_t Length() const { return m_ * n_; }
  

  /* 
   * \brief Resizes 2D array  
   * \warning existing data not preserved
   * \todo maybe another implementation that preserves existing data by
   */
  void Resize(const size_t &m, const size_t &n) {
    m_ = m;
    n_ = n;
    data_.resize(m_*n_);
  }

  void Resize(const size_t &m, const size_t &n, const T &t) {
    data_.clear();
    m_ = m;
    n_ = n;
    data_.resize(m_ *n_, t);
  }

  /// \brief Set all values in array to given value
  void SetValue(const T &t) {
    for(size_t i = 0; i < data_.size(); ++i) {
      data_[i] = t;
    }
  }

  /**
   * \brief Assign / access Array2D values (row-major order)
   *
   * If "my_data" is an object of type Array2D, then its array values are
   * my_data(i,j) = value_ij.
   */
  T& operator()(size_t i, const size_t j) { return data_[j + n_ * i]; }
  const T& operator()(size_t i,size_t j) const { return data_[j + n_ * i]; }

  private:
    size_t m_; ///> Number of elements in 1st-dimension
    size_t n_; ///> Number of elements in 2nd-dimension

    /**
     * \brief Data in array with size = m_ * n_
     *
     * The data is stored in row-major order (2nd index is fast-varying).
     */
    vector<T> data_;
};

#endif // ARRAY2D_H_
