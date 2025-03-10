### OpenSBLI CFD solver HPC performance benchmark (2024) ###
### Authors: Dr. D.J. Lusher (Japan Aerospace Exploration Agency, JAXA), Prof. N.D. Sandham (University of Southampton). Additional credit to Prof S.P. Jammy (SRM University) for an earlier version of this application on which this benchmark is based. ###

### This application is provided as a computational benchmark for performance and scaling tests for the code on different types of HPC architectures. ###

### If this application is used to present performance data, please cite the following OpenSBLI code references: ###
[1] D.J. Lusher, S.P. Jammy, N.D. Sandham. OpenSBLI: Automated code-generation for heterogeneous computing architectures applied to compressible fluid dynamics on structured grids. Computer Physics Communications 267 (2021).
[2] D.J.Lusher, A. Sansica, N.D. Sandham, J. Meng, B. Siklosi, A.Hashimoto. OpenSBLI v3.0: High-Fidelity Multi-Block Transonic Aerofoil CFD Simulations using Domain Specific Languages on GPUs. Computer Physics Communications (2024).

### Parallel execution is provided courtesy of the Oxford Parallel Structured (OPS) domain specific language for structured mesh computations (https://github.com/OP-DSL/OPS). Please reference the following for further details: ###
[3] I. Z. Reguly, G. R. Mudalige and M. B. Giles, Loop Tiling in Large-Scale Stencil Codes at Run-Time with OPS, in IEEE Transactions on Parallel and Distributed Systems, vol. 29, no. 4, pp. 873-886, 1 April 2018, doi: 10.1109/TPDS.2017.2778161.

### The full-3D unsteady Navier-Stokes equations are solved for the supersonic Taylor-Green vortex test case introduced in: ###
[4] D.J. Lusher, N.D. Sandham. Assessment of low-dissipative shock-capturing schemes for the compressible Taylor–Green vortex. AIAA-Journal, 59 (2), 533-545 (2021).
[5] J.B. Chapelier, D.J. Lusher, et al. Comparison of high-order numerical methodologies for the simulation of the supersonic Taylor–Green vortex flow. Physics of Fluids 36 (5) (2024).

### Instructions ###

### The test case is provided with a choice of two different numerical methods: ###

1. TENO6: Direct application of a 6th-order Targeted Essentially Non-Oscillatory shock-capturing scheme applied to the convective terms.

2. Central4-WENO5: A 4th-order central difference scheme in the Kennedy-Gruber-Pirozolli (KGP) cubic-split form, with 5th-order Weighted Essentially Non-Oscillatory shock-capturing applied as a filter-step method.

In both cases, 4th-order central differences are used to compute diffusive and heat-conduction terms. A 3rd-order low-storage Runge-Kutta scheme is used for time-advancement. Please refer to the main OpenSBLI references [1,2] for further details.


### Configuring the input file: ###
Each directory contains a file called 'input', which controls the options for the performance benchmark. Control of the number of CPU/GPU resources, MPI ranks, and OpenMP threads should be specified in the job submission script depending on the machine being used.

The input file has the following options in order:
1. Simulation type: va, ss, ws (see below)
2. Number of grid points per direction.
3. Number of time-step iterations to perform.
4. Number of MPI processes (weak-scaling only).
5. OPS detailed kernel-by-kernel timing output: 0 or 1.
6. HDF5 file I/O timing.

For example, for a weak-scaling test with N=64^3 points per process, with 4 processes and 100 iterations, no kernel timing output nor HDF5 I/O timing:

ws 64 100 4 0 False

### Validation simulation (va): ###
To check the correct execution of the code on a system prior to performing performance and scaling tests, a validation case can be run to check numerical output values. Please select the "va" option in the input file. If the simulation has executed correctly, the following numerical outputs should be obtained.

Iteration, Time, KE, dilatation_dissipation, enstrophy_dissipation, rhom
500, 8.46250e-01, 2.65363e+05, 1.16991e+05, 1.79666e+06, 2.09715e+06

The validation simulation always performs 500 iterations on a fixed N=128^3 mesh.

### Strong scaling test (ss): ###
Select the "ss" option to perform strong scaling tests on a fixed number of mesh points. Adjust the number of CPUs/GPUs in the job submission script.

### Weak scaling test (ss): ###
Select the "ws" option to perform weak scaling tests on a fixed number of mesh points per process. Adjust the number of CPUs/GPUs in the job submission script.



























