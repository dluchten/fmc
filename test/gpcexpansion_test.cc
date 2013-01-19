#include "gpcexpansion.h"
#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <cstdio>
#include <iostream>

// Cannot change these parameters without updating test functions
const int order1 = 20, num_dims1 = 1;
const int order2 = 26, num_dims2 = 2;

class GPCExpansionTest : public testing::Test {
protected:
  virtual void SetUp() {
    myexp1_ = new GPCExpansion("Legendre", order1, num_dims1);
    myexp2_ = new GPCExpansion("Legendre", order2, num_dims2);
  }
  GPCExpansion *myexp1_;
  GPCExpansion *myexp2_;
};

/*
 * Test 1D expansion
 *
 */
TEST_F(GPCExpansionTest, OneDimensionalExpansionIsCorrect) {
  double kTol = 1.e-14;

  // Compute Expansion for sin(x) on [a,b]
  double a = 0., b = 1.1 * acos(-1.);
  double xmid = (a + b) / 2.;
  double dx2 = (b - a) / 2.;
  int num_nodes = myexp1_->num_nodes();
  Array1D<double> nodes(num_nodes);
  nodes = myexp1_->oned_nodes();
  for (int i = 0; i != num_nodes; ++i) {
    double x = xmid + dx2 * nodes(i);
    double value  = sin(x);
    myexp1_->set_nodal_value(i, value);
  }
  myexp1_->CompCoefficients();

  // Eval expansion and check error
  int num_points = 100;
  double dx = (b - a) / (num_points - 1);
  double x = a;
  double error_max = 1.e-20;
  Array1D<double> point(1);
  double value;
  for (int i = 0; i != num_points; ++i) {
    point(0) = (x - a) / dx2 - 1.;
    myexp1_->EvalExpansion(point, value);
    double exact_value = sin(x);
    double error = exact_value - value;
    error_max = max(error, error_max);
    x += dx;
  }
  ASSERT_NEAR(0., error_max, kTol);
}


/*
 * Test 2D expansion
 *
 */
TEST_F(GPCExpansionTest, TwoDimensionalExpansionIsCorrect) {
  double kTol = 1.e-14;

  // Compute Expansion for sin(x) * sin(y) on [a,b] x [c,d]
  double a = 0., b = 1. * acos(-1.);
  double c = 0., d = 1.33 * acos(-1);
  double xmid = (a + b) / 2.;
  double dx2 = (b - a) / 2.;
  double ymid = (c + d) / 2.;
  double dy2 = (d - c) / 2.;
  int num_oned_nodes = myexp2_->num_oned_nodes();
  Array1D<double> oned_nodes(num_oned_nodes);
  oned_nodes = myexp2_->oned_nodes();
  for (int i = 0; i != num_oned_nodes; ++i) {
    for (int j = 0; j != num_oned_nodes; ++j) {
    double x = xmid + dx2 * oned_nodes(j);
    double y = ymid + dy2 * oned_nodes(i);
    double value  = sin(x) * sin(y);
    int index = j + i * num_oned_nodes;
    myexp2_->set_nodal_value(index, value);
    }
  }
  myexp2_->CompCoefficients();

  // Eval expansion and check error
  int num_points = 100;
  double dx = (b - a) / (num_points - 1);
  double dy = (d - c) / (num_points - 1);
  double x = a, y = c;
  double error_max = 1.e-20;
  Array1D<double> point(2);
  double value;
  for (int i = 0; i != num_points; ++i) {
    x = a;
    for (int j = 0; j != num_points; ++j) {
    point(0) = (x - a) / dx2 - 1.;
    point(1) = (y - c) / dy2 - 1.;
    myexp2_->EvalExpansion(point, value);
    double exact_value = sin(x) * sin(y);
    double error = exact_value - value;

    // printf("(x, y) = % 6.4f % 6.4f  ", x, y);

    error_max = max(error, error_max);
    x += dx;
    }

    // printf("\n");

    y += dy;
  }
  ASSERT_NEAR(0., error_max, kTol);
}
