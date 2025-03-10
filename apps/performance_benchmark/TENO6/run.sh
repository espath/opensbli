# Regular execution
mpirun -n 1 opensbli_mpi_cuda
mpirun -n 2 opensbli_mpi_cuda
mpirun -n 4 opensbli_mpi_cuda
# GPU-direct
#mpirun -n 1 opensbli_mpi_cuda
#mpirun -n 2 opensbli_mpi_cuda -gpudirect
#mpirun -n 4 opensbli_mpi_cuda -gpudirect
# Tiled
mpirun -n 1 opensbli_mpi_cuda_tiled OPS_TILING OPS_TILING_MAXDEPTH=9 OPS_TILESIZE_X=10000 OPS_TILESIZE_Y=10000 OPS_TILESIZE_Z=10000
mpirun -n 2 opensbli_mpi_cuda_tiled OPS_TILING OPS_TILING_MAXDEPTH=9 OPS_TILESIZE_X=10000 OPS_TILESIZE_Y=10000 OPS_TILESIZE_Z=10000
mpirun -n 4 opensbli_mpi_cuda_tiled OPS_TILING OPS_TILING_MAXDEPTH=9 OPS_TILESIZE_X=10000 OPS_TILESIZE_Y=10000 OPS_TILESIZE_Z=10000
