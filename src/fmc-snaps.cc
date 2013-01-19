#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "orthopolybasis.h"
#include "legendrebasis.h"
#include "gpcexpansion.h"
#include "grid2d.h"
#include "array1d.h"
#include "velocityfield_snap_2d.h"
#include "doublegyre.h"
#include "integrator.h"
#include "heun.h"
#include "runge-kutta.h"
#include "flowmap2d.h"
#include "stopwatch.h"
#include "typedefs.h"

int main(int argc, char *argv[]) {

  int vis_setting = 0;
  int mgrid, ngrid;
  if ((argc < 2) || (argc > 4)) {
    printf ("USAGE: %s [index_file] [vis_options]\n", argv[0]);
    printf ("  [vis_options] is (blank) for default\n");
    printf ("                   \"snap\" for use resolution of snaps\n");
    printf ("                   nx ny to set # of pixels manually\n");
    exit(EXIT_FAILURE);
  }
  if (argc == 3) {
    if (strcmp(argv[2], "snap") == 0) {
      vis_setting = 1;
    } else {
      printf ("USAGE: %s [index_file] [vis_options]\n", argv[0]);
      printf ("  [vis_options] is (blank) for default\n");
      printf ("                   \"snap\" for use resolution of snaps\n");
      printf ("                   nx ny to set # of pixels manually\n");
      exit(EXIT_FAILURE);
    }
  }
  if (argc == 4) {
    vis_setting = 2;
    mgrid = atoi(argv[2]);
    ngrid = atoi(argv[3]);
  }

  // Parse the index file
  // FILE LIST MUST BE AT LEAST 2 FILES LONG
  // Format of file is a list of filenames followed by times
  // e.g. file.0 0.00
  //      file.1 0.01
  //      file.3 0.02
  //      ...    ...
  // NOTE: Coordinates taken from only the first of these files
  //       Grid assumed to be static (and regular) throughout files
  //       Whitespace ignored between filenames and times
  // TODO: Add errors for reading files that are bad
  std::ifstream ifs;
  std::vector<std::string> flist;
  VecDoub1D time_vec;
  ifs.open(argv[1]);
  assert(ifs.good());

  std::string temp_str;
  double temp_doub;
  while( ifs >> temp_str >> temp_doub) {
    flist.push_back(temp_str);
    time_vec.push_back(temp_doub);
  }

  // Start stopwatch
  Stopwatch stopwatch;
  
  // Init snapshot velocity field, and set bounds to those in files
  VelocityFieldSnap2D *vf = new VelocityFieldSnap2D(time_vec, flist);
  double xmin = vf->xlo();
  double xmax = vf->xhi();
  double ymin = vf->ylo();
  double ymax = vf->yhi();

  // Guess an appropriate dt and num_steps, given the input time vector
  double dt = (time_vec.back() - time_vec[0]) / (time_vec.size()-1);
  int num_steps = time_vec.size()-1;
  
  // TODO grid viz mgrid and ngrid interpreted from user?

  // Grid for viz.
  if (vis_setting == 0) { // default
    mgrid = 400;
    ngrid = 200;
  } else if (vis_setting == 1) { // snap
    mgrid = vf->nx();
    ngrid = vf->ny();
  } else if (vis_setting == 2) { // Xpix Ypix
    ;
  } else {
    std::cout << "Bad vis setting" << std::endl;
    exit(EXIT_FAILURE);
  }
    
  Grid2D grid(mgrid, ngrid, xmin, xmax, ymin, ymax);
  Grid2D vals = grid;

  int order = 20;
  FlowMap2D flowmap(order, xmin, xmax, ymin, ymax);
  Grid2D nodes; // quadrature nodes
  Grid2D nodes_dt; // quadrature nodes after dt
  flowmap.GetNodes(nodes);
  flowmap.GetNodes(nodes_dt);

  // Init Timestepper
  Integrator *integrator = new RungeKutta4(dt, *vf);

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



  // Deallocate
  delete integrator;
  delete vf;

  stopwatch.toc();
  return 0;
}
