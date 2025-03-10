make clean
rm *.h5
make -j12 opensbli_mpi_cuda
make -j12 opensbli_mpi_cuda_tiled
