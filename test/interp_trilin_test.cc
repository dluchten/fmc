#include "interp_trilin.h"
#include <gtest/gtest.h>

class InterpTrilinTest : public testing::Test {
protected:
  virtual void SetUp() {
    // Testing 2 functions over the same 3x4x5 spatial grid
    for (int i=0; i<3; ++i) x.push_back(0.);
    for (int i=0; i<4; ++i) y.push_back(0.);
    for (int i=0; i<5; ++i) z.push_back(0.);
    f0xyz = new VecDoub3D(3, VecDoub2D(4, VecDoub1D(5, 0.)));
    f1xyz = new VecDoub3D(3, VecDoub2D(4, VecDoub1D(5, 0.)));

    // Set values for the coordinates and for the 2 functions
    for (int i=0; i<3; ++i) x[i] =  4. + 1.*i;
    for (int i=0; i<4; ++i) y[i] =  9. + 1.*i;
    for (int i=0; i<5; ++i) z[i] = -1. + 1.*i;
    for (int i=0; i<3; ++i) {
      for (int j=0; j<4; ++j) {
	for (int k=0; k<5; ++k) {
	  (*f0xyz)[i][j][k] = 7.*x[i]*x[i] + 11.*y[j] + 13.*z[k];
	  (*f1xyz)[i][j][k] = 2.*x[i] - 21.*x[i]*y[j] +  4.*z[k]*z[k];
	}
      }
    }

    // Create trilinear interpolator, initially interpolating over f0xyz
    mytlint = new InterpTrilin(x, y, z, f0xyz);

    // Space to hold coordinate triples temporarily
    xyz = new double[3];
  }
  VecDoub1D x, y, z;         // spatial coordinate vectors
  VecDoub3D *f0xyz, *f1xyz;  // 3D arrays (functions) over x,y,z coord vectors
  InterpTrilin *mytlint;     // trilinear interpolator
  double *xyz;               // space to hold temporary coord triples
  double fxyz;               // space to hold interpolated function values
};

// Interpolate non-linear function over the given grid
// Points selected are halfway between grid points
// Therefore, the interpolated value is the average of
// the function values on the grid points
TEST_F(InterpTrilinTest, InterpolationIsCorrect) {
  xyz[0] =  4.5;
  xyz[1] = 11.5;
  xyz[2] = -0.5;
  // f0xyz(4.5, 11.5, -0.5) = 7*4.5*4.5 + 11*11.5 + 13*(-0.5) = 261.75
  // Using x values of 4 and 5
  // Using y values of 11 and 12
  // Using z values of -1 and 0
  // f0xyz(4, 11, -1) = 220
  // f0xyz(4, 11,  0) = 233
  // f0xyz(4, 12, -1) = 231
  // f0xyz(4, 12,  0) = 244
  // f0xyz(5, 11, -1) = 283
  // f0xyz(5, 11,  0) = 296
  // f0xyz(5, 12, -1) = 294
  // f0xyz(5, 12,  0) = 307
  // Interpolated value = 263.5
  fxyz = mytlint->Interpolate(xyz);
  ASSERT_DOUBLE_EQ(263.5, fxyz);
}

// Interpolators must be able to change the field over which they interpolate
// We use this functionality to interpolate vector fields
// (multiple scalar fields over the same grid)
TEST_F(InterpTrilinTest, AbleToChangePointer) {
  xyz[0] =  4.5;
  xyz[1] = 11.5;
  xyz[2] = -0.5;
  // f1xyz(4.5, 11.5, -0.5) = 2*4.5 - 21*4.5*11.5 + 4*(-0.5)*(-0.5) = -1076.75
  // f1xyz(4, 11, -1) = - 912
  // f1xyz(4, 11,  0) = - 916
  // f1xyz(4, 12, -1) = - 996
  // f1xyz(4, 12,  0) = -1000
  // f1xyz(5, 11, -1) = -1141
  // f1xyz(5, 11,  0) = -1145
  // f1xyz(5, 12, -1) = -1246
  // f1xyz(5, 12,  0) = -1250
  // Interpolated value = -1075.75
  mytlint->ChangeFieldPtr(f1xyz);
  fxyz = mytlint->Interpolate(xyz);
  ASSERT_DOUBLE_EQ(-1075.75, fxyz);
}
