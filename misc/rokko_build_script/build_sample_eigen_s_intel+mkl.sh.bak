#!/bin/bash -x

rm -R ../../build/rokko/*

cmake ~/development/rokko/ -DEIGEN_S_LIB="/home/sakashita/eigen_s/libEigen_s.a -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" -DEIGEN_S_INC=$HOME/eigen_s -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_Fortran_COMPILER="mpif90" -DCMAKE_BUILD_TYPE=Debug

make diagonalize_eigen_s VERBOSE=1

#make diagonalize_time_frank_eigen_s VERBOSE=1
