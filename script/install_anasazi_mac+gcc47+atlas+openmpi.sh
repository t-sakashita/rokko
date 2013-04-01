#!/bin/bash

PREFIX="$1"
test -z "$PREFIX" && PREFIX=/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf trilinos-11.0.3-Source
tar jxf $HOME/source/trilinos-11.0.3-Source.tar.bz2 

mkdir trilinos-11.0.3-build
cd trilinos-11.0.3-build

cmake -DTPL_ENABLE_MPI=ON -DCMAKE_Fortran_COMPILER="openmpif90" -DCMAKE_CXX_COMPILER="openmpicxx" -DCMAKE_C_COMPILER="openmpicc" -DTPL_BLAS_LIBRARIES="-L/opt/local -lcblas -lf77blas -latlas" -DTPL_LAPACK_LIBRARIES="-L/opt/local -llapack" -DCMAKE_INSTALL_PREFIX="$PREFIX" -DTrilinos_ENABLE_Anasazi=ON $HOME/build/trilinos-11.0.3-Source

make -j4
make install
