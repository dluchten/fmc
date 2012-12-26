#include "legendrebasis.h"
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include <cstdio>
#include <iostream>

class LegendreBasisTest : public testing::Test {
protected:
  virtual void SetUp() {
    myleg4 = new LegendreBasis(4);
    myleg5 = new LegendreBasis(5);
  }
  GPCBasis *myleg4;
  GPCBasis *myleg5;
};

/*
 * Test Quadrature Rule - InitNodesWeights()
 *
 * Take exact values from mathworld, note the symmetry
 * http://mathworld.wolfram.com/LobattoQuadrature.html
 */
TEST_F(LegendreBasisTest, QuadratureRuleIsCorrect) {
  // TEST 1: quadrature rule for 5 nodes
  int num_nodes = 4 + 1; // interpolating polynomial: (order + 1) nodes
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
  ASSERT_EQ( num_nodes, myleg4->num_nodes());
  for (int i = 0; i < num_nodes; ++i) {
    ASSERT_DOUBLE_EQ(exact_nodes.at(i), myleg4->nodes(i));
    ASSERT_DOUBLE_EQ(exact_weights.at(i), myleg4->weights(i));
  }

  // TEST 2: 6 quadrature nodes
  num_nodes = 5 + 1;  // interpolating polynomial: (order + 1) nodes
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
  ASSERT_EQ(num_nodes, myleg5->num_nodes());
  for (int i = 0; i < num_nodes; ++i) {
    ASSERT_DOUBLE_EQ(exact_nodes.at(i), myleg5->nodes(i));
    ASSERT_DOUBLE_EQ(exact_weights.at(i), myleg5->weights(i));
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
  int num_nodes = myleg4->num_nodes();
  int order = myleg4->order();
  ASSERT_DOUBLE_EQ(4, order);
  std::vector<double> p(order + 1); // store evalution of basis functions

  double z; // store node i in z
  // Check all nodes i
  for (int i = 0; i < num_nodes; ++i) {
    z = myleg4->nodes(i);
    p[0] = 1.;
    p[1] = z;
    p[2] = 1. / 2. * (3. * z * z - 1.);
    p[3] = 1. / 2. * (5. * pow(z, 3) - 3. * z);
    p[4] = 1. / 8. * (35. * pow(z, 4) - 30. * z * z + 3.);
    // Check all degrees j
    for (int j = 0; j < order + 1; ++j ) {
      ASSERT_DOUBLE_EQ(p.at(j), myleg4->basisfunctionsatnodes(i,j));
    }
  }
}
