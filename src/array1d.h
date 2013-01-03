#ifndef ARRAY1D_H_
#define ARRAY1D_H_

/**
 * \file array1d.h
 * \brief Class for storage of 1D array (vector) of specified type T
 * \author D. M. Luchtenburg
 *
 * Implementation of 1D array to store data.
 */
#include <vector>
#include <cstddef> 
using namespace std;

/**
 * \brief Stores data in 1D array (vector) of specified type T
 */
template <typename T>
class Array1D {
 public:
  /// \brief Default constructor (no allocation)
  Array1D() : n_(0) {}

  /// \brief Constructor that allocates memory
  Array1D(const size_t &n) : n_(n) { data_.resize(n_); }

  /** \brief Constructor that allocates memory, and initializes data to constant
   * value
   */
  Array1D(const size_t &n, const T &t) : n_(n) { data_.resize(n_, t); }

  /// \brief Destructor
  ~Array1D() { data_.clear(); }

  /* /// \brief Function to clear the memory */
  /* void Clear() { n_ = 0; data_.clear(); } */

  /// \brief Returns length
  size_t Length() const { return n_; }

  /// \brief Resizes the array
  void Resize(const size_t &n) { n_ = n; data_.resize(n_); }

  /// \brief Resizes the array and sets values
  void Resize(const size_t &n, const T &t) {
    data_.clear();
    n_ = n;
    data_.resize(n_, t);
  }

  /// \brief Set all values in array to given value 
  void SetValue(const T &t) {
    for(size_t i = 0; i < data_.size(); ++i) {
      data_[i] = t;
    }
  }

  /**
   * \brief Assign / access Array1D values
   *
   * If "my_data" is an object of type Array1D, then its array values are
   * my_data(i) = value
   */
  T& operator()(size_t i) { return data_[i]; } ///< assign
  const T& operator()(size_t i) const { return data_[i]; } ///< read

    // Neat concept: returns reference to data_[i] and 
    // this can directly assigned or read

  /// \brief Add element t to the end of the vector
  void PushBack(const T &t) { n_ += 1; data_.push_back(t); }

 private:
  size_t n_; 
  vector<T> data_;
};

#endif // ARRAY1D_H_
