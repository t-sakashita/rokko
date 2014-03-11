#!/bin/bash -x

rm -R ../../build/rokko/*

cmake $WORK/development/rokko/ \
-DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort \
-DCMAKE_INSTALL_PREFIX="$WORK/lib/rokko" \
-DMPI_C_INCLUDE_PATH="/usr/include" -DMPI_CXX_INCLUDE_PATH="/usr/include" \
-DMPI_C_LIBRARIES="-lmpi" -DMPI_CXX_LIBRARIES="-lmpi++;-lmpi" -DMPI_Fortran_LIBRARIES="-lmpi" \
-DSCALAPACK_LIB="-lmkl_scalapack_lp64 -lmkl_blacs_sgimpt_lp64 -lifcore -mkl=parallel" \
-DCMAKE_BUILD_TYPE=Debug \
-DBOOST_INCLUDE_DIR="/opt/nano/alps/boost_1_52_0/"


make diagonalize_scalapack VERBOSE=1

#make diagonalize_time_frank_eigen_s VERBOSE=1
