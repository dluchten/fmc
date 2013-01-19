#ifndef _FLOWMAP2D_
#define _FLOWMAP2D_

/** 
 * \file flowmap2d.h
 * \brief flow map for 2D.
 * \author D. M. Luchtenburg
 *
 * This class describes the implementation of the 2D flow mapthe Legendre basis functions.
 */
#include "gpcexpansion.h"
#include "array1d.h"
#include "array2d.h"
#include "grid2d.h"

class FlowMap2D {
public:
  FlowMap2D(int order, double xmin, double xmax, double ymin, double ymax);
  ~FlowMap2D();
  /// \brief Get nodes
  void GetNodes(Grid2D &nodes);
  /// \brief Set nodal values
  void SetNodalValues(const Grid2D &nodalvals);
  /// \brief Compute expansion coefficients
  void CompCoefficients();

  /// \brief Evaluates flow map for grid vals.
  void Eval(Grid2D &vals) const;

  /// \brief Compute Jacobian
  void CompJac(const Grid2D &val);

  /// \brief Initialize FTLE field
  void InitFTLE(const int mgrid, const int ngrid);
  
  /// \brief Initialize FTLE field
  void CompFTLE(const double &t, const Grid2D &val);

  /// \brief Initialize FTLE field
  void DumpFTLE(const char* filename) const;

  private:
  /// \brief Initialize Jacobian
  void InitJac(const int m, const int n);

  /// \brief Compute eigenvalues of 2x2 matrix
  void EigenValues(const Array2D<double> &mat, Array1D<double> &eigvals) const;

  /// \brief Compute max. singular value of 2x2 matrix
  double MaxSingularValue(const Array2D<double> &mat) const;

    /// \brief Convert coordinates to coordinates on standard element [-1, 1]
  void CoordToStandardEl(const double &x, const double &y,
                         double &xse, double &yse) const;
  /// \brief Convert coordinates on standard element [-1, 1] to coordinates
  void StandardElToCoord(const double &xse, const double &yse,
                         double &x, double &y) const;

  /// \brief Scaling values for conversion coordinates <-> standard element [-1, 1]
  double dx_, dx_half_, x_mid_;
  double dy_, dy_half_, y_mid_;
  
  /// \brief Order of flow map expansion.
  int order_;
  /// \brief Interval for x-values: [xmin_, xmax_].
  double xmin_, xmax_;
  /// \brief Interval for y-values: [ymin_, ymax_].
  double ymin_, ymax_;
  /// \brief flow map expansion for x coordinate.
  GPCExpansion *flowmap_xdir_;
  /// \brief flow map expansion for x coordinate.
  GPCExpansion *flowmap_ydir_;

  /// \brief Dimensions of viz grid array
  int mgrid_, ngrid_;

  // \brief viz grid
  Grid2D grid;

  /// \brief store Jacobian matrices
  Array2D<double> D11_old_;
  Array2D<double> D12_old_;
  Array2D<double> D21_old_;
  Array2D<double> D22_old_;
  Array2D<double> D11_new_;
  Array2D<double> D12_new_;
  Array2D<double> D21_new_;
  Array2D<double> D22_new_;

  /// \brief store FTLE field
  Array2D<double> ftle_;
};

#endif  //  _FLOWMAP2D_

