#include "legendrebasis.h"
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include <cstdio>
#include <iostream>

// Cannot change these parameters without updating test functions
// (analytical values are hardcoded)
const int order1 = 4, num_nodes1 = 5;
const int order2 = 3, num_nodes2 = 6;

class LegendreBasisTest : public testing::Test {
protected:
  virtual void SetUp() {
    myleg1_ = new LegendreBasis(order1, num_nodes1); // order, num_nodes
    myleg2_ = new LegendreBasis(order2, num_nodes2);
  }
  OrthoPolyBasis *myleg1_;
  OrthoPolyBasis *myleg2_;
};

/*
 * Test Quadrature Rule - InitNodesWeights()
 *
 * Take exact values from mathworld, note the symmetry
 * http://mathworld.wolfram.com/LobattoQuadrature.html
 */
TEST_F(LegendreBasisTest, QuadratureRuleIsCorrect) {
  // TEST 1: quadrature rule for 5 nodes
  int num_nodes = myleg1_->num_nodes();
  ASSERT_EQ(num_nodes1, num_nodes);
  std::vector<double> exact_nodes(num_nodes);
  std::vector<double> exact_weights(num_nodes);
  // Quadrature nodes
  exact_nodes[0] = -1.; 
  exact_nodes[num_nodes - 1] = -exact_nodes[0];
  exact_nodes[1] = -1. / 7. * sqrt(21.); 
  exact_nodes[num_nodes - 2] = -exact_nodes[1];
  exact_nodes[2] = 0.;
  // Quadrature weights
  exact_weights[0] = 1. / 10.; 
  exact_weights[num_nodes - 1] = exact_weights[0];
  exact_weights[1] = 49. / 90.;
  exact_weights[num_nodes - 2] = exact_weights[1];
  exact_weights[2] = 32. / 45.;
  // Compare exact values with computed ones
  for (int i = 0; i < num_nodes; ++i) {
    ASSERT_DOUBLE_EQ(exact_nodes.at(i), myleg1_->nodes(i));
    ASSERT_DOUBLE_EQ(exact_weights.at(i), myleg1_->weights(i));
  }

  // TEST 2: 6 quadrature nodes
  num_nodes = myleg2_->num_nodes();
  ASSERT_EQ(num_nodes2, num_nodes);
  exact_nodes.resize(num_nodes); 
  exact_weights.resize(num_nodes);
  // Quadrature nodes
  exact_nodes[0] = -1.; 
  exact_nodes[num_nodes - 1] = -exact_nodes[0];
  exact_nodes[1] = -sqrt(1. / 21. * (7 + 2. * sqrt(7.))); 
  exact_nodes[num_nodes - 2] = -exact_nodes[1];
  exact_nodes[2] = -sqrt(1. / 21. * (7 - 2. * sqrt(7.))); 
  exact_nodes[num_nodes - 3] = -exact_nodes[2];
  // Quadrature weights
  exact_weights[0] = 1. / 15.; 
  exact_weights[num_nodes - 1] = exact_weights[0];
  exact_weights[1] = 1. / 30. * (14. - sqrt(7.));
  exact_weights[num_nodes - 2] = exact_weights[1];
  exact_weights[2] = 1. / 30. * (14. + sqrt(7.));
  exact_weights[num_nodes - 3] = exact_weights[2];
  // Compare exact values with computed ones
  for (int i = 0; i < num_nodes; ++i) {
    ASSERT_DOUBLE_EQ(exact_nodes.at(i), myleg2_->nodes(i));
    ASSERT_DOUBLE_EQ(exact_weights.at(i), myleg2_->weights(i));
  }
}

/*
 * Test Evaluate Basis Functions at nodes - EvalBasisFunctionsAtNodes()
 *
 * Take exact values from mathworld, note the symmetry
 * http://mathworld.wolfram.com/LobattoQuadrature.html
 * http://mathworld.wolfram.com/LegendrePolynomial.html
 */
TEST_F(LegendreBasisTest, EvalBasisFunctionsAtNodesIsCorrect) {
  // TEST 1
  int num_nodes = myleg1_->num_nodes();
  int order = myleg1_->order();
  ASSERT_EQ(order1, order);
  ASSERT_EQ(num_nodes1, num_nodes);
  std::vector<double> p(order + 1); // basis functions (const. - order)
  double z; // store node i in z
  // Check all nodes i
  for (int i = 0; i < num_nodes; ++i) {
    z = myleg1_->nodes(i);
    p[0] = 1.;
    p[1] = z;
    p[2] = 1. / 2. * (3. * z * z - 1.);
    p[3] = 1. / 2. * (5. * pow(z, 3) - 3. * z);
    p[4] = 1. / 8. * (35. * pow(z, 4) - 30. * z * z + 3.);
    // Check all degrees j
    for (int j = 0; j < order + 1; ++j ) {
      ASSERT_DOUBLE_EQ(p.at(j), myleg1_->basisfunctionsatnodes(i,j));
    }
  }
  // TEST 2
  num_nodes = myleg2_->num_nodes();
  order = myleg2_->order();
  ASSERT_EQ(order2, order);
  ASSERT_EQ(num_nodes2, num_nodes);
  p.resize(order + 1);
  // Check all nodes i
  for (int i = 0; i < num_nodes; ++i) {
    z = myleg2_->nodes(i);
    p[0] = 1.;
    p[1] = z;
    p[2] = 1. / 2. * (3. * z * z - 1.);
    p[3] = 1. / 2. * (5. * pow(z, 3) - 3. * z);
    // Check all degrees j
    for (int j = 0; j < order + 1; ++j ) {
      ASSERT_DOUBLE_EQ(p.at(j), myleg2_->basisfunctionsatnodes(i,j));
    }
  }
}


TEST_F(LegendreBasisTest, TestOrthogonalityBasisFunctions) {
  // TEST 1
  int num_nodes = myleg2_->num_nodes();
  int order = myleg2_->order();
  ASSERT_EQ(order2, order);
  ASSERT_EQ(num_nodes2, num_nodes);
  const double kTol = 1.e-15;
  double p; // store inner product in p

  // Check all basis functions i
  for (int i = 0; i < order + 1; ++i) {
  // Check all basis functions j
    for (int j = i; j < order + 1; ++j) {
      // discrete inner product (quad. rule) 
      p = 0.;
      for (int k = 0; k < num_nodes; ++k) {
        p += myleg2_->basisfunctionsatnodes(k, i) *
            myleg2_->basisfunctionsatnodes(k, j) *
            myleg2_->weights(k);
      }
      if (i == j) { 
        ASSERT_DOUBLE_EQ(2. / (2. * i + 1.), p);
      } else {
        ASSERT_NEAR(0., p, kTol);
      }
    }
  }  
}
