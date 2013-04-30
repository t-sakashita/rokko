#!/bin/bash -x

rm -R ../../build/rokko/*

cmake ~/development/rokko/  -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_Fortran_COMPILER="mpif90" -DCMAKE_BUILD_TYPE=Debug -DELEMENTAL_DIR="/opt/nano/rokko"

make diagonalize_elemental VERBOSE=1

make diagonalize_time_frank_elemental VERBOSE=1
