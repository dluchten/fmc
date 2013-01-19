#include "velocityfield_snap_2d.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <assert.h>

VelocityFieldSnap2D::VelocityFieldSnap2D(const VecDoub1D &t,
					 const VecDoub1D &x, const VecDoub1D &y)
  : VelocityFieldSnap(2) {

  // Check if t, x and y are monotonically increasing
  for (size_t i=1; i<x.size(); ++i) assert(x[i] > x[i-1]);
  for (size_t i=1; i<y.size(); ++i) assert(y[i] > y[i-1]);
  for (size_t i=1; i<t.size(); ++i) assert(t[i] > t[i-1]);

  t_ = new VecDoub1D(t);
  x_ = new VecDoub1D(x);
  y_ = new VecDoub1D(y);
  velx_ = new VecDoub3D(t.size(), VecDoub2D(x.size(), VecDoub1D(y.size(), 0.)));
  vely_ = new VecDoub3D(t.size(), VecDoub2D(x.size(), VecDoub1D(y.size(), 0.)));
  interp_ = new InterpTrilin(*t_, *x_, *y_, NULL);
  coords_ = new double[3];
}

VelocityFieldSnap2D::VelocityFieldSnap2D(const VecDoub1D &t,
					 const std::vector<std::string> &flist)
  : VelocityFieldSnap(2) {

  // Check that the number of files in the file list matches the number of times
  if (flist.size() != t.size()) {
    std::cout << "# of input files must match # of times" << std::endl;
    exit(1);
  }

  std::ifstream ifs;

  // Get grid information from the first file in the list
  // First must parse the header
  ifs.open(flist[0].c_str());
  assert(ifs.good());
  VecDoub1D x, y;
  std::string temp;

  std::getline(ifs, temp, '=');
  if (temp != "VARIABLES ") {
    std::cout << "Could not find VARIABLES in input file header" << std::endl;
    exit(1);
  }
  std::getline(ifs, temp);
  if (temp != " \"X\", \"Y\", \"U\", \"V\"") {
    std::cout << "No \"X\", \"Y\", \"U\", \"V\" in header" << std::endl;
    exit(1);
  }
  std::getline(ifs, temp, '=');
  if (temp != "ZONE I") {
    std::cout << "Could not find ZONE I in input file header" << std::endl;
    exit(1);
  }
  std::getline(ifs, temp, ',');
  x.resize(atoi(temp.c_str()));
  std::getline(ifs, temp, '=');
  if (temp != " J") {
    std::cout << "Could not find a J= in the input file header" << std::endl;
    exit(1);
  }
  std::getline(ifs, temp, ',');
  y.resize(atoi(temp.c_str()));
  std::getline(ifs, temp);
  if (temp != " DATAPACKING=BLOCK") {
    std::cout << "Must be DATAPACKING=BLOCK" << std::endl;
    exit(1);
  }

  // Now get the grid coordinates ... first x, then y
  for (size_t i=0; i<x.size(); ++i) {
    ifs >> x[i];
  }
  std::getline(ifs, temp);
  double tempdoub;
  for (size_t i=0; i<y.size(); ++i) {
    ifs >> y[i];
    for (int j=0; j<(x.size()-1); ++j) {
      ifs >> tempdoub;
    }
  }
  ifs.close();

  // Check that all vectors are monotonically increasing
  for (size_t i=1; i<x.size(); ++i) assert(x[i] > x[i-1]);
  for (size_t i=1; i<y.size(); ++i) assert(y[i] > y[i-1]);
  for (size_t i=1; i<t.size(); ++i) assert(t[i] > t[i-1]);

  // Construct velocity field using read data
  t_ = new VecDoub1D(t);
  x_ = new VecDoub1D(x);
  y_ = new VecDoub1D(y);
  velx_ = new VecDoub3D(t.size(), VecDoub2D(x.size(), VecDoub1D(y.size(), 0.)));
  vely_ = new VecDoub3D(t.size(), VecDoub2D(x.size(), VecDoub1D(y.size(), 0.)));
  interp_ = new InterpTrilin(*t_, *x_, *y_, NULL);
  coords_ = new double[3];

  // Add the velocities from all of the files
  // Assume that the grid is the same for all files!
  int indices[2];
  for (size_t i=0; i<flist.size(); ++i) {
    ifs.open (flist[i].c_str());
    assert(ifs.good());

    // Skip header and grid coordinates
    std::getline(ifs, temp);
    std::getline(ifs, temp);
    std::getline(ifs, temp);
    std::getline(ifs, temp);

    // Read in the x velocities
    for (size_t j=0; j<y.size(); ++j) {
      indices[1] = j;
      for (size_t k=0; k<x.size(); ++k) {
	indices[0] = k;
	ifs >> tempdoub;
	Set(i, &indices[0], 0, tempdoub);
      }
    }

    // Read in the y velocities
    for (size_t j=0; j<y.size(); ++j) {
      indices[1] = j;
      for (size_t k=0; k<x.size(); ++k) {
	indices[0] = i;
	ifs >> tempdoub;
	Set(i, &indices[0], 0, tempdoub);
      }
    }
    ifs.close();
  }
}

VelocityFieldSnap2D::~VelocityFieldSnap2D() {
  delete velx_;
  delete vely_;
  delete t_;
  delete x_;
  delete y_;
  delete interp_;
  delete [] coords_;
}

void VelocityFieldSnap2D::Get(double t, const double *x, double *vf) const {
  coords_[0] = t;
  coords_[1] = x[0];
  coords_[2] = x[1];

  // Interpolate x-velocity
  interp_->ChangeFieldPtr(velx_);
  vf[0] = interp_->Interpolate(&coords_[0]);

  // Interpolate y-velocity
  interp_->ChangeFieldPtr(vely_);
  vf[1] = interp_->Interpolate(&coords_[0]);
}

int VelocityFieldSnap2D::Set(const int t, const int *x,
			     const int comp, const double vcomp) {

  // Check bounds
  int num_pts[2];
  num_pts[0] = x_->size();
  num_pts[1] = y_->size();
  const int snap_capacity = (*velx_).size();

  for (int i=0; i<2; ++i) {
    if ((x[i] < 0) || (x[i] > (num_pts[i]-1))) {
	std::cout << "Out of bounds for velocity insertion" << std::endl;
	exit(1);
    }
  }
  if ((t < 0) || (t > (snap_capacity-1))) {
    std::cout << "Out of bounds time ind for velocity insertion" << std::endl;
    exit(1);
  }

  if (comp == 0) {
    (*velx_)[t][x[0]][x[1]] = vcomp;
  } else if (comp == 1) {
    (*vely_)[t][x[0]][x[1]] = vcomp;
  } else {    
    std::cout << "Bad component for velocity insertion" << std::endl;
    exit(1);
  }
  return 0;
}
