#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "grid2d.h"
#include "array1d.h"
#include "doublegyre.h"
#include "integrator.h"
#include "heun.h"
#include "runge-kutta.h"
#include "flowmap2d.h"
#include "stopwatch.h"

using std::string;

/**
 * \brief Main program
 *
 * Only 2D supported
 */

int main(int argc, char *argv[]) {

  // Handle input string
  if (argc != 7) {
    printf("USAGE: %s <integrator> <timestep> <numsteps> <order> <mgrid> <ngrid>\n", argv[0]);
    printf(" with the following choices:\n");
    printf(" <integrator>: heun, rk4\n");
    exit(EXIT_FAILURE);
  }
  const string integrator_string = argv[1]; // type of integrator for time-stepping
  const double dt = atof(argv[2]);          // time step
  const int num_steps = atoi(argv[3]);      // number of time steps
  const int order = atoi(argv[4]);          // order of expansion
  const int mgrid = atoi(argv[5]);          // grid size for visualization
  const int ngrid = atoi(argv[6]);

  // Start stopwatch
  Stopwatch stopwatch;
  
  // Init double gyre velocity field on [0, 2] x [0, 1]
  VelocityField *model = new DoubleGyre();
  double xmin = 0., xmax = 2.;
  double ymin = 0., ymax = 1.;

  // Grid for visualization
  Grid2D grid(mgrid, ngrid, xmin, xmax, ymin, ymax);
  Grid2D vals = grid;

  // Define short-time flow map 
  FlowMap2D flowmap(order, xmin, xmax, ymin, ymax);
  Grid2D nodes;    // quadrature nodes
  Grid2D nodes_dt; // quadrature nodes after dt
  // Get quadrature nodes
  flowmap.GetNodes(nodes); flowmap.GetNodes(nodes_dt);

  // Initialize Timestepper
  Integrator *integrator;
  if (integrator_string.compare("heun") == 0) {
    integrator = new Heun(dt, *model);
  } else if (integrator_string.compare("rk4") == 0) {
    integrator = new RungeKutta4(dt, *model);
  } else {
    cout << "Invalid choice for <integrator>" << endl;
    exit(EXIT_FAILURE);
  }

  // Initialize FTLE
  flowmap.InitFTLE(mgrid, ngrid);
  
  // Compute flow map and FTLE
  // using composition
  double t = 0.; 
  for (int i = 0; i < num_steps; ++i) {
    integrator->Step(t, nodes, nodes_dt);
    t = (i + 1) * dt;
    // Short-time flow map approximation
    flowmap.SetNodalValues(nodes_dt);
    flowmap.CompCoefficients();
    // Evaluate flow map on visualization grid
    flowmap.Eval(vals);
    // Compute FTLE on visualization grid
    flowmap.CompFTLE(t, vals);
    // Write FTLE field
    char filename[256];
    sprintf(filename,"ftle_%03d.dat", i+1);
    flowmap.DumpFTLE(filename);
  }

  // Deallocation
  delete integrator;
  delete model;

  // Time calculation
  stopwatch.toc();
  return 0;
}
