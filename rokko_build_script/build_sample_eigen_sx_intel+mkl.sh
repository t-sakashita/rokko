#!/bin/bash -x

rm -R ../../build/rokko/*

cmake ~/development/rokko/ -DEIGEN_SX_LIB="/home/sakashita/eigen_sx/libEigen_sx.a -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" -DEIGEN_SX_INC=$HOME/eigen_sx -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_Fortran_COMPILER="mpif90" -DCMAKE_BUILD_TYPE=Debug

make diagonalize_eigen_sx VERBOSE=1

##make diagonalize_time_frank_eigen_sx VERBOSE=1
