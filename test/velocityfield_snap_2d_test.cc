#include "velocityfield_snap_2d.h"
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>

class VelocityFieldSnap2DTest : public testing::Test {
protected:
  virtual void SetUp() {
    nt1 = 12;
    nx1 = 100;
    ny1 = 200;

    // Resize all the vectors
    t1.resize(nt1);
    x1.resize(nx1);
    y1.resize(ny1);
    vx1.resize(nt1);
    vy1.resize(nt1);
    for (int i=0; i<nt1; ++i) vx1[i].resize(nx1);
    for (int i=0; i<nt1; ++i) vy1[i].resize(nx1);
    for (int i=0; i<nt1; ++i) {
      for (int j=0; j<nx1; ++j) {
	vx1[i][j].resize(ny1);
      }
    }
    for (int i=0; i<nt1; ++i) {
      for (int j=0; j<nx1; ++j) {
	vy1[i][j].resize(ny1);
      }
    }

    // Set values for the coordinates
    for (int i=0; i<nt1; ++i) t1[i] = 0.1*i - 0.5;
    for (int i=0; i<nx1; ++i) x1[i] = i;
    for (int i=0; i<ny1; ++i) y1[i] = ny1 - i;

    // Set arbitrary velocity field using arbitrary functions of t,x,y
    for (int i=0; i<nt1; ++i) {
      for (int j=0; j<nx1; ++j) {
	for (int k=0; k<ny1; ++k) {
	  vx1[i][j][k] = t1[i] * (x1[j]+y1[k]);
	  vy1[i][j][k] = t1[i]*x1[j] - y1[k];
	}
      }
    }

    vfsnap1 = new VelocityFieldSnap2D(t1, x1, y1);

    // Objects needed for vfsnap2 (construct later in test)
    nt2 = 4;
    t2.resize(nt2);
    t2[0] = 0.0;
    t2[1] = 0.1;
    t2[2] = 0.2;
    t2[3] = 0.3;
    flist2.resize(nt2);

    // List of small files to read through (Tecplot output)
    flist2[0] = "velocityfield_snap_2d_test.0.dat";
    flist2[1] = "velocityfield_snap_2d_test.1.dat";
    flist2[2] = "velocityfield_snap_2d_test.2.dat";
    flist2[3] = "velocityfield_snap_2d_test.3.dat";
  }
  VelocityFieldSnap2D *vfsnap1;     // velocityfield using data from Set
  VelocityFieldSnap2D *vfsnap2;     // velocityfield using data from files
  int nx1, ny1;                     // size grid in x,y for vfsnap1
  int nt1, nt2;                     // size of time grid for vfsnap{1,2}
  VecDoub1D x1, y1, t1, t2;         // time and space vectors
  VecDoub3D vx1, vy1;               // gridded velocities for vfsnap1
  std::vector<std::string> flist2;  // list of input files for vfsnap2
};

// The vfsnap1 in memory gets a dimension upon construction
// Since we have 2 spatial dimensions, it should return 2
TEST_F(VelocityFieldSnap2DTest, ReturnsDimension) {
  int dim = vfsnap1->dimen();
  ASSERT_EQ(2, dim);
}

// Here we use the Set function to "fill" vfsnap1 with velocity data
// We are looping through the arrays vx1 and vy1 to fill with these values
// After the arrays are full, ask for the function values at the grid points
// These should be exactly the same as the input data, as long as the
// function values do not vary tremendously between the grid points
TEST_F(VelocityFieldSnap2DTest, FillsSnaps) {

  // Load the arrays into the VelocityFieldSnap2D class
  int coord_ind[2];
  for (int i=0; i<nt1; ++i) {
    for (int j=0; j<nx1; ++j) {
      coord_ind[0] = j;
      for (int k=0; k<ny1; ++k) {
	coord_ind[1] = k;
	vfsnap1->Set(i, &coord_ind[0], 0, vx1[i][j][k]);
      }
    }
    for (int j=0; j<nx1; ++j) {
      coord_ind[0] = j;
      for (int k=0; k<ny1; ++k) {
	coord_ind[1] = k;
	vfsnap1->Set(i, &coord_ind[0], 1, vy1[i][j][k]);
      }
    }
  }
  
  // Verify that the contents have been copied
  int vx_orig, vy_orig;
  int vx_copy, vy_copy;
  double coord_dbl[2], v_copy[2];
  for (int i=0; i<nt1; ++i) {
    for (int j=0; j<nx1; ++j) {
      coord_dbl[0] = x1[j];
      for (int k=0; k<ny1; ++k) {
	coord_dbl[1] = y1[k];
	vfsnap1->Get(t1[i], &coord_dbl[0], &v_copy[0]);
	vx_orig = vx1[i][j][k];
	vy_orig = vy1[i][j][k];
	vx_copy = v_copy[0];
	vy_copy = v_copy[1];
	ASSERT_DOUBLE_EQ(vx_orig, vx_copy);
	ASSERT_DOUBLE_EQ(vy_orig, vy_copy);
      }
    }
  }
}

// vfsnap2 is constructed by reading a series of files
// Test passes if we are able to construct it
TEST_F(VelocityFieldSnap2DTest, ConstructsFromSeriesOfTecplotFiles) {
    vfsnap2 = new VelocityFieldSnap2D(t2, flist2);
}

// Use constructed vfsnap2 from last test
// Read through the same files used to construct vfsnap2,
// and compare with velocity values on the grid points
// These should match
TEST_F(VelocityFieldSnap2DTest, InterpolatesCorrectlyFromTecplotFiles) {

  double vel_first[2];
  double vel_interpolated[2];
  double pos_first[2];
  std::string temp;

  std::ifstream ifs;
  ifs.open(flist2[0].c_str());

  // Read beyond header
  std::getline(ifs, temp);
  std::getline(ifs, temp);

  // Get first x and first y positions
  ifs >> pos_first[0];
  std::getline(ifs, temp);
  ifs >> pos_first[1];
  std::getline(ifs, temp);
  ifs >> vel_first[0];
  std::getline(ifs, temp);
  ifs >> vel_first[1];
  ifs.close();

  // See if the returned velocity matches that from the file
  vfsnap2->Get(t2[0], &pos_first[0], &vel_interpolated[0]);
  for (int i=0; i<2; ++i) {
    ASSERT_DOUBLE_EQ(vel_first[i], vel_interpolated[i]);
  }
}
