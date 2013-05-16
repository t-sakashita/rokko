#!/bin/bash -x

rm -R ../../build/rokko/*

cmake ~/development/rokko/ \
-DSCALAPACK_LIB="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" \
-DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_Fortran_COMPILER="mpif90" -DCMAKE_BUILD_TYPE=Debug

make diagonalize_eigen_s VERBOSE=1

#make diagonalize_time_frank_eigen_s VERBOSE=1
