#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "orthopolybasis.h"
#include "gpcexpansion.h"

/**
 * \brief Main program
 *
 * \todo Maybe change data structures using boost
 */

int main(int argc, char *argv[]) {

  int order = 2, num_dims = 2;
  GPCExpansion *myexp = new GPCExpansion("Legendre", order, num_dims);


  // // Test test basis


  // GPCBasis *mybasis = new LegendreBasis(3);


  // int num_nodes = mybasis->num_nodes();
  // int order = mybasis->order();
  // double z; 

  // std::cout << "num_nodes = " << num_nodes << "\n";
  // std::cout << "order = " << order << "\n";



  // z = mybasis->basisfunctionsatnodes(num_nodes-1, order);
  // cout << "z = " << z;

  
  // // Handle input string
  // ParmParser parser(argc, argv);

  // // // Select 1) ODE (provide model), or 2) snapshots (provide file with list)
  // // string velocity_field_selector = parser.GetString(
  // //     "velocity_field_selector", "Specification of velocity field (ODE, snapshots)", "ODE");
  // // string velocity_field_def = parser.GetString(
  // //     "velocity_field_def", "ODE model (doublegyre), or file with list of snapshots", "doublegyre");

  // // ONLY SNAPSHOTS IMPLEMENTED FOR NOW
  // string fname_snap_list = parser.GetString(
  //     "fname_snap_list", "file with list of snapshots", "snapshots.dat");
  // // t0, te, dt
  // double tstart = parser.GetDouble(
  //     "tstart", "start time for computation", 0.);
  // double tend = parser.GetDouble(
  //     "tend", "end time for computation", 10.);
  // double dt = parser.GetDouble(
  //     "dt", "time step for integration", .1);
  // // xmin, xmax, ymin, ymax
  // double xmin = parser.GetDouble(
  //     "xmin", "lower bound for x of grid domain", 0.);
  // double xmax = parser.GetDouble(
  //     "xmax", "upper bound for x of grid domain", 2.);
  // double ymin = parser.GetDouble(
  //     "ymin", "lower bound for y of grid domain", 0.);
  // double ymax = parser.GetDouble(
  //     "ymax", "upper bound for y of grid domain", 1.);
  // // nx, ny, type flow map basis
  // int nx = parser.GetInt(
  //     "nx", "number of points in x-direction", 200);
  // int ny = parser.GetInt(
  //     "ny", "number of points in y-direction", 100);
  // string flow_map_type = parser.GetString(
  //     "flow_map_type", "flow map basis type (nneigbor,leg)","nneigbor" );
  // // interpolator type
  // string interpolator_type = parser.GetString(
  //     "interpolator_type", "interpolator type (nneigbor,spline)","nneigbor");
  // // time stepper
  // string time_stepper_type = parser.GetString(
  //     "time_stepper_type", "Timestepper type (Heun)","Heun");
  // // outputdir
  // string outputdir = parser.GetString(
  //     "outputdir", "Output directory", "./out/");

  // // Set snapshot grid 
  // SnapGrid snapgrid(fnamesnap_list);

  // // Set particle grid (non necessarily same as snapshot grid)
  // ParticleGrid grid(nx,ny,xmin,xmax,ymin,ymax,flow_map_type);

  // // Init spatial interpolator for velocity field at non snapshot grid points
  // Interpolator *interp = InitInterpolator(snapgrid,interpolator_type);

  // // Init velocity field 
  // string velocity_field_type = "snapshots";
  // VelocityField *velocity_field = InitVelocityField(
  //     velocity_field_type, fnamesnap_list,interp);

  // // Init Timestepper
  // Integrator* integrator = InitIntegrator(time_stepper_type, grid, 
  //                                         *velocity_field);
  // // Init short-time flow maps
  // Flowmap flowmaps[nsteps] = new FlowMap(flow_map_type,xmin,xmax,ymin,ymax,nx,ny);

  // // Integrate for spec. # steps
  // for(int i = 0; i < nsteps; ++i) {
  //   integrator->Step(t,grid.get_particles(),gridout);
  //   t = (i + 1) * dt;
  //   // Compute short time flow map
  //   flowmaps[i]->set_flowmap(gridout);
  //   } 
  //   timer.print();

  // delete mybasis;
  return 0;
}
