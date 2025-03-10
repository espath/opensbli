#include <stdlib.h> 
#include <string.h> 
#include <math.h> 
#include "constants.h"
#define OPS_3D
#define OPS_API 2
#include "ops_seq.h"
#include "opensbliblock00_kernels.h"
#include "io.h"
#include "reductions.h"

int main(int argc, char **argv) 
{
// Initializing OPS 
char simulation_type[2], HDF5op[5];
int number_of_processes, grid_points, number_of_iters, OPSdiags;
FILE *inputfile;
inputfile = fopen("input", "r");
if (inputfile) {
fscanf(inputfile, "%s %d %d %d %d %s", simulation_type, &grid_points, &number_of_iters, &number_of_processes, &OPSdiags, HDF5op);
}
// Check whether OPS kernel-based timing output is required
if (OPSdiags == 1){
ops_init(argc,argv,2);
}
else{
ops_init(argc,argv,1);
}
ops_printf("\n-------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
ops_printf("This is the OpenSBLI V3.0 (2024) CFD benchmark application developed for the UK HPC benchmark suite. \nAuthors: Dr. David J. Lusher (JAXA), Prof. Satya P. Jammy (SRM University), and Prof. Neil Sandham (University of Southampton).\n\n");
ops_printf("Further information on the code and numerical methods: (Computer Physics Communications journal (2021, 2024)):\n \nD.J. Lusher et al. OpenSBLI: Automated code-generation for heterogeneous computing architectures applied to compressible fluid dynamics on structured grids.\nD.J. Lusher et al. OpenSBLI v3.0: High-Fidelity Multi-Block Transonic Aerofoil CFD Simulations using Domain Specific Languages on GPUs.");
ops_printf("\n-------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
ops_printf("Inputs: Simulation: %s, Grid points: %d, Iterations: %d, Processes (weak-scaling mode only): %d, OPS kernel output: %d, HDF5 output: %s", simulation_type, grid_points, number_of_iters, number_of_processes, OPSdiags, HDF5op);
ops_printf("\n-------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n");

if (strcmp(simulation_type,"va") == 0){
ops_printf("Validation simulation selected, defaults to N=128^3 mesh.\n");
block0np0 = 128;
block0np1 = 128;
block0np2 = 128;
Lx0 = 2.0*M_PI;
Lx1 = 2.0*M_PI;
Lx2 = 2.0*M_PI;
niter = 500;
print_iter = 100;
dt = 0.00084625;
dt = dt*2.0;
post_process = 1;
ops_printf("The number of grid points used per direction is N = %d\nThe total number of grid points is N = %e\n",block0np0, (double) block0np0*block0np1*block0np2);
}
else if (strcmp(simulation_type,"ss") == 0){
ops_printf("Strong scaling simulation selected.\n");
block0np0 = grid_points;
block0np1 = grid_points;
block0np2 = grid_points;
post_process = 0;
Lx0 = 2.0*M_PI;
Lx1 = 2.0*M_PI;
Lx2 = 2.0*M_PI;
niter = number_of_iters;
print_iter = niter/10;

dt = 0.00084625*256/grid_points;
ops_printf("The number of grid points used per direction is %d\nThe total number of grid points is %e\n",block0np0, (double) block0np0*block0np1*block0np2);
}
else if (strcmp(simulation_type,"ws") == 0){
ops_printf("Weak scaling simulation selected on %d MPI processes using N = %d^3 grid points per process \n", number_of_processes, grid_points);
block0np0 = grid_points*number_of_processes;
block0np1 = grid_points*number_of_processes;
block0np2 = grid_points*number_of_processes;

post_process = 0;
Lx0 = 2.0*M_PI*number_of_processes;
Lx1 = 2.0*M_PI*number_of_processes;
Lx2 = 2.0*M_PI*number_of_processes;
niter = number_of_iters;
print_iter = niter/10;

dt = 0.00084625*256/grid_points;
ops_printf("The total number of grid points is: N = %e\n", (double) block0np0*block0np1*block0np2);
}
else{
ops_printf("The simulation name is invalid \nThe simulation name provided is %s\nit should be one of the following\nva, ss or ws\n",simulation_type);
ops_printf("Exiting the simulation\n");
exit(0);
}
// Set restart to 1 to restart the simulation from HDF5 file
restart = 0;
Delta0block0 = Lx0/block0np0;
Delta1block0 = Lx1/block0np1;
Delta2block0 = Lx2/block0np2;
eps = 1.00000000000000e-40;
TENO_CT = 1e-6;
double rkB[] = {(1.0/3.0), (15.0/16.0), (8.0/15.0)};
double rkA[] = {0, (-5.0/9.0), (-153.0/128.0)};
write_output_file = 10000;
HDF5_timing = 1;
Minf = 1.25;
gama = 1.4;
Re = 1600.0;
Pr = 0.71;
shock_filter_control = 1.00000000000000;
gamma_m1 = -1 + gama;
inv2Delta0block0 = 1.0/(Delta0block0*Delta0block0);
inv2Delta1block0 = 1.0/(Delta1block0*Delta1block0);
inv2Delta2block0 = 1.0/(Delta2block0*Delta2block0);
inv2Minf = 1.0/(Minf*Minf);
invDelta0block0 = 1.0/(Delta0block0);
invDelta1block0 = 1.0/(Delta1block0);
invDelta2block0 = 1.0/(Delta2block0);
invPr = 1.0/(Pr);
invRe = 1.0/(Re);
inv_gamma_m1 = 1.0/((-1 + gama));
invgamma_m1 = 1.0/(gamma_m1);
ops_decl_const("Delta0block0" , 1, "double", &Delta0block0);
ops_decl_const("Delta1block0" , 1, "double", &Delta1block0);
ops_decl_const("Delta2block0" , 1, "double", &Delta2block0);
ops_decl_const("HDF5_timing" , 1, "int", &HDF5_timing);
ops_decl_const("Minf" , 1, "double", &Minf);
ops_decl_const("Pr" , 1, "double", &Pr);
ops_decl_const("Re" , 1, "double", &Re);
ops_decl_const("TENO_CT" , 1, "double", &TENO_CT);
ops_decl_const("block0np0" , 1, "int", &block0np0);
ops_decl_const("block0np1" , 1, "int", &block0np1);
ops_decl_const("block0np2" , 1, "int", &block0np2);
ops_decl_const("dt" , 1, "double", &dt);
ops_decl_const("eps" , 1, "double", &eps);
ops_decl_const("gama" , 1, "double", &gama);
ops_decl_const("gamma_m1" , 1, "double", &gamma_m1);
ops_decl_const("inv2Delta0block0" , 1, "double", &inv2Delta0block0);
ops_decl_const("inv2Delta1block0" , 1, "double", &inv2Delta1block0);
ops_decl_const("inv2Delta2block0" , 1, "double", &inv2Delta2block0);
ops_decl_const("inv2Minf" , 1, "double", &inv2Minf);
ops_decl_const("invDelta0block0" , 1, "double", &invDelta0block0);
ops_decl_const("invDelta1block0" , 1, "double", &invDelta1block0);
ops_decl_const("invDelta2block0" , 1, "double", &invDelta2block0);
ops_decl_const("invPr" , 1, "double", &invPr);
ops_decl_const("invRe" , 1, "double", &invRe);
ops_decl_const("inv_gamma_m1" , 1, "double", &inv_gamma_m1);
ops_decl_const("invgamma_m1" , 1, "double", &invgamma_m1);
ops_decl_const("niter" , 1, "int", &niter);
ops_decl_const("shock_filter_control" , 1, "double", &shock_filter_control);
ops_decl_const("simulation_time" , 1, "double", &simulation_time);
ops_decl_const("start_iter" , 1, "int", &start_iter);
ops_decl_const("write_output_file" , 1, "int", &write_output_file);
// Define and Declare OPS Block
ops_block opensbliblock00 = ops_decl_block(3, "opensbliblock00");
#include "defdec_data_set.h"
// Define and declare stencils
#include "stencils.h"
// Define and declare OPS reduction handles
double dilatation_dissipation_B0_out = 0.0;
ops_reduction dilatation_dissipation_B0 = ops_decl_reduction_handle(sizeof(double), "double", "reduction_dilatation_dissipation_B0");
double KE_B0_out = 0.0;
ops_reduction KE_B0 = ops_decl_reduction_handle(sizeof(double), "double", "reduction_KE_B0");
double enstrophy_dissipation_B0_out = 0.0;
ops_reduction enstrophy_dissipation_B0 = ops_decl_reduction_handle(sizeof(double), "double", "reduction_enstrophy_dissipation_B0");
double rhom_B0_out = 0.0;
ops_reduction rhom_B0 = ops_decl_reduction_handle(sizeof(double), "double", "reduction_rhom_B0");
#include "bc_exchanges.h"
// Init OPS partition
double partition_start0, elapsed_partition_start0, partition_end0, elapsed_partition_end0;
ops_timers(&partition_start0, &elapsed_partition_start0);
ops_partition("");
ops_timers(&partition_end0, &elapsed_partition_end0);
ops_printf("-----------------------------------------\n MPI partition and reading input file time: %lf\n -----------------------------------------\n", elapsed_partition_end0-elapsed_partition_start0);
// Constants from HDF5 restart file
if (restart == 1){
ops_get_const_hdf5("simulation_time", 1, "double", (char*)&simulation_time, "restart.h5");
ops_get_const_hdf5("iter", 1, "int", (char*)&start_iter, "restart.h5");
}
else {
simulation_time = 0.0;
start_iter = 0;
}
tstart = simulation_time;

if (restart == 0){
int iteration_range_43_block0[] = {-5, block0np0 + 5, -5, block0np1 + 5, -5, block0np2 + 5};
ops_par_loop(opensbliblock00Kernel043, "Grid_based_initialisation0", opensbliblock00, 3, iteration_range_43_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(x0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(x1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(x2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_idx());
}

// Initialize loop timers
double cpu_start0, elapsed_start0, cpu_end0, elapsed_end0;
ops_timers(&cpu_start0, &elapsed_start0);
double inner_start, elapsed_inner_start;
double inner_end, elapsed_inner_end;
ops_timers(&inner_start, &elapsed_inner_start);

ops_halo_transfer(periodicBC_direction0_side0_37_block0);
ops_halo_transfer(periodicBC_direction0_side1_38_block0);
ops_halo_transfer(periodicBC_direction1_side0_39_block0);
ops_halo_transfer(periodicBC_direction1_side1_40_block0);
ops_halo_transfer(periodicBC_direction2_side0_41_block0);
ops_halo_transfer(periodicBC_direction2_side1_42_block0);

for(iter=start_iter; iter<=start_iter+niter - 1; iter++)
{
simulation_time = tstart + dt*((iter - start_iter)+1);
if(fmod(iter+1, print_iter) == 0){
        ops_timers(&inner_end, &elapsed_inner_end);
        ops_printf("Iteration: %d. Time-step: %.3e. Simulation time: %.5f. Time/iteration: %lf.\n", iter+1, dt, simulation_time, (elapsed_inner_end - elapsed_inner_start)/print_iter);
        ops_NaNcheck(rho_B0);
        ops_timers(&inner_start, &elapsed_inner_start);
}

for(stage=0; stage<=2; stage++)
{
int iteration_range_3_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel003, "CRu0", opensbliblock00, 3, iteration_range_3_block0,
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_11_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel011, "CRu2", opensbliblock00, 3, iteration_range_11_block0,
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_12_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel012, "CRu1", opensbliblock00, 3, iteration_range_12_block0,
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_5_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel005, "CRp", opensbliblock00, 3, iteration_range_5_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_6_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel006, "CRa", opensbliblock00, 3, iteration_range_6_block0,
ops_arg_dat(p_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(a_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_15_block0[] = {-2, block0np0 + 2, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel015, "CRT_B0", opensbliblock00, 3, iteration_range_15_block0,
ops_arg_dat(p_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(T_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_14_block0[] = {-2, block0np0 + 2, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel014, "CRmu_B0", opensbliblock00, 3, iteration_range_14_block0,
ops_arg_dat(T_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(mu_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_0_block0[] = {-1, block0np0 + 1, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel000, "LFTeno_reconstruction_0_direction", opensbliblock00, 3, iteration_range_0_block0,
ops_arg_dat(a_B0, 1, stencil_0_01_00_00_5, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_23_00_00_13, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_23_00_00_13, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_23_00_00_13, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_23_00_00_13, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_23_00_00_13, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_23_00_00_13, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_23_00_00_13, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_01_00_00_5, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_01_00_00_5, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk3_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk4_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_1_block0[] = {0, block0np0, -1, block0np1 + 1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel001, "LFTeno_reconstruction_1_direction", opensbliblock00, 3, iteration_range_1_block0,
ops_arg_dat(a_B0, 1, stencil_0_00_01_00_5, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_00_23_00_13, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_00_23_00_13, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_23_00_13, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_23_00_13, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_23_00_13, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_23_00_13, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_00_01_00_5, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_00_23_00_13, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_00_01_00_5, "double", OPS_READ),
ops_arg_dat(wk5_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk6_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk7_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk8_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk9_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_2_block0[] = {0, block0np0, 0, block0np1, -1, block0np2 + 1};
ops_par_loop(opensbliblock00Kernel002, "LFTeno_reconstruction_2_direction", opensbliblock00, 3, iteration_range_2_block0,
ops_arg_dat(a_B0, 1, stencil_0_00_00_01_5, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_00_00_23_13, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_23_13, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_23_13, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_23_13, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_23_13, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_23_13, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_00_00_01_5, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_00_00_01_5, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_00_00_23_13, "double", OPS_READ),
ops_arg_dat(wk10_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk11_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk12_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk13_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(wk14_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_13_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel013, "LFTeno Residual", opensbliblock00, 3, iteration_range_13_block0,
ops_arg_dat(wk0_B0, 1, stencil_0_10_00_00_5, "double", OPS_READ),
ops_arg_dat(wk10_B0, 1, stencil_0_00_00_10_5, "double", OPS_READ),
ops_arg_dat(wk11_B0, 1, stencil_0_00_00_10_5, "double", OPS_READ),
ops_arg_dat(wk12_B0, 1, stencil_0_00_00_10_5, "double", OPS_READ),
ops_arg_dat(wk13_B0, 1, stencil_0_00_00_10_5, "double", OPS_READ),
ops_arg_dat(wk14_B0, 1, stencil_0_00_00_10_5, "double", OPS_READ),
ops_arg_dat(wk1_B0, 1, stencil_0_10_00_00_5, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_10_00_00_5, "double", OPS_READ),
ops_arg_dat(wk3_B0, 1, stencil_0_10_00_00_5, "double", OPS_READ),
ops_arg_dat(wk4_B0, 1, stencil_0_10_00_00_5, "double", OPS_READ),
ops_arg_dat(wk5_B0, 1, stencil_0_00_10_00_5, "double", OPS_READ),
ops_arg_dat(wk6_B0, 1, stencil_0_00_10_00_5, "double", OPS_READ),
ops_arg_dat(wk7_B0, 1, stencil_0_00_10_00_5, "double", OPS_READ),
ops_arg_dat(wk8_B0, 1, stencil_0_00_10_00_5, "double", OPS_READ),
ops_arg_dat(wk9_B0, 1, stencil_0_00_10_00_5, "double", OPS_READ),
ops_arg_dat(Residual0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual3_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE),
ops_arg_dat(Residual4_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_16_block0[] = {0, block0np0, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel016, "Derivative evaluation CD u0_B0 x0 ", opensbliblock00, 3, iteration_range_16_block0,
ops_arg_dat(u0_B0, 1, stencil_0_22_00_00_8, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_18_block0[] = {0, block0np0, -2, block0np1 + 2, 0, block0np2};
ops_par_loop(opensbliblock00Kernel018, "Derivative evaluation CD u1_B0 x0 ", opensbliblock00, 3, iteration_range_18_block0,
ops_arg_dat(u1_B0, 1, stencil_0_22_00_00_8, "double", OPS_READ),
ops_arg_dat(wk1_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_20_block0[] = {0, block0np0, 0, block0np1, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel020, "Derivative evaluation CD u2_B0 x0 ", opensbliblock00, 3, iteration_range_20_block0,
ops_arg_dat(u2_B0, 1, stencil_0_22_00_00_8, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_22_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel022, "Derivative evaluation CD u0_B0 x1 ", opensbliblock00, 3, iteration_range_22_block0,
ops_arg_dat(u0_B0, 1, stencil_0_00_22_00_8, "double", OPS_READ),
ops_arg_dat(wk3_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_23_block0[] = {0, block0np0, 0, block0np1, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel023, "Derivative evaluation CD u1_B0 x1 ", opensbliblock00, 3, iteration_range_23_block0,
ops_arg_dat(u1_B0, 1, stencil_0_00_22_00_8, "double", OPS_READ),
ops_arg_dat(wk4_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_24_block0[] = {0, block0np0, 0, block0np1, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel024, "Derivative evaluation CD u2_B0 x1 ", opensbliblock00, 3, iteration_range_24_block0,
ops_arg_dat(u2_B0, 1, stencil_0_00_22_00_8, "double", OPS_READ),
ops_arg_dat(wk5_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_25_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel025, "Derivative evaluation CD u0_B0 x2 ", opensbliblock00, 3, iteration_range_25_block0,
ops_arg_dat(u0_B0, 1, stencil_0_00_00_22_8, "double", OPS_READ),
ops_arg_dat(wk6_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_26_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel026, "Derivative evaluation CD u1_B0 x2 ", opensbliblock00, 3, iteration_range_26_block0,
ops_arg_dat(u1_B0, 1, stencil_0_00_00_22_8, "double", OPS_READ),
ops_arg_dat(wk7_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_27_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel027, "Derivative evaluation CD u2_B0 x2 ", opensbliblock00, 3, iteration_range_27_block0,
ops_arg_dat(u2_B0, 1, stencil_0_00_00_22_8, "double", OPS_READ),
ops_arg_dat(wk8_B0, 1, stencil_0_00_00_00_3, "double", OPS_WRITE));

int iteration_range_35_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel035, "Viscous terms", opensbliblock00, 3, iteration_range_35_block0,
ops_arg_dat(T_B0, 1, stencil_0_22_22_22_27, "double", OPS_READ),
ops_arg_dat(mu_B0, 1, stencil_0_22_22_22_27, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_22_22_22_27, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_22_22_22_27, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_22_22_22_27, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_00_22_22_19, "double", OPS_READ),
ops_arg_dat(wk1_B0, 1, stencil_0_00_22_00_11, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_00_00_22_11, "double", OPS_READ),
ops_arg_dat(wk3_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk4_B0, 1, stencil_0_00_00_22_11, "double", OPS_READ),
ops_arg_dat(wk5_B0, 1, stencil_0_00_00_22_11, "double", OPS_READ),
ops_arg_dat(wk6_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk7_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(wk8_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(Residual1_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(Residual2_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(Residual3_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(Residual4_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW));

int iteration_range_49_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel049, "Temporal solution advancement", opensbliblock00, 3, iteration_range_49_block0,
ops_arg_dat(Residual0_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(Residual1_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(Residual2_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(Residual3_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(Residual4_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhoE_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rho_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou0_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou1_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_dat(rhou2_RKold_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW),
ops_arg_gbl(&rkA[stage], 1, "double", OPS_READ),
ops_arg_gbl(&rkB[stage], 1, "double", OPS_READ));

ops_halo_transfer(periodicBC_direction0_side0_37_block0);
ops_halo_transfer(periodicBC_direction0_side1_38_block0);
ops_halo_transfer(periodicBC_direction1_side0_39_block0);
ops_halo_transfer(periodicBC_direction1_side1_40_block0);
ops_halo_transfer(periodicBC_direction2_side0_41_block0);
ops_halo_transfer(periodicBC_direction2_side1_42_block0);
}
// End of the time loop
}

if (post_process == 1){
ops_printf("\nValidation simulation: The output below should match that in the README file.\n");
int iteration_range_48_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel048, "Taylor-Green vortex post-processing", opensbliblock00, 3, iteration_range_48_block0,
ops_arg_dat(mu_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00_00_00_3, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_22_22_22_27, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_22_22_22_27, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_22_22_22_27, "double", OPS_READ),
ops_arg_reduce(KE_B0, 1, "double", OPS_INC),
ops_arg_reduce(dilatation_dissipation_B0, 1, "double", OPS_INC),
ops_arg_reduce(enstrophy_dissipation_B0, 1, "double", OPS_INC),
ops_arg_reduce(rhom_B0, 1, "double", OPS_INC),
ops_arg_dat(divV_B0, 1, stencil_0_00_00_00_3, "double", OPS_RW));

ops_reduction_result(enstrophy_dissipation_B0, &enstrophy_dissipation_B0_out);
ops_reduction_result(dilatation_dissipation_B0, &dilatation_dissipation_B0_out);
ops_reduction_result(KE_B0, &KE_B0_out);
ops_reduction_result(rhom_B0, &rhom_B0_out);
// Print output
ops_printf("Iteration, Time, KE, dilatation_dissipation, enstrophy_dissipation, rhom\n");
// Write the output values
ops_printf("%d, %.5e, %.5e, %.5e, %.5e, %.5e\n", iter, simulation_time, KE_B0_out, dilatation_dissipation_B0_out, enstrophy_dissipation_B0_out, rhom_B0_out);
}

ops_timers(&cpu_end0, &elapsed_end0);
ops_printf("\nTimings are:\n");
ops_printf("-----------------------------------------\n");
ops_printf("Total Wall time %lf\n",elapsed_end0-elapsed_start0);

if (strcmp(HDF5op,"True") == 0){
HDF5_IO_Write_0_opensbliblock00(opensbliblock00, rho_B0, rhou0_B0, rhou1_B0, rhou2_B0, rhoE_B0, x0_B0, x1_B0, x2_B0, HDF5_timing);
}
else{
ops_printf("Not performing file I/O as write HDF5 is set to False\n");
}
if (OPSdiags == 1){
ops_printf("OPS kernel performance output:\n");
ops_timing_output(std::cout);
}
ops_exit();
//Main program end 
}
