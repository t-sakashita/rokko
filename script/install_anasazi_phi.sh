#!/bin/bash

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf trilinos-11.4.1-Source
tar jxf $HOME/source/trilinos-11.4.1-Source.tar.bz2

mkdir trilinos-11.4.1-build
rm -rf trilinos-11.4.1-build/*
cd trilinos-11.4.1-build

cmake \
-D TPL_ENABLE_MPI:BOOL=ON -D MPI_USE_COMPILER_WRAPPERS:BOOL=ON \
-D MPI_Fortran_COMPILER:FILEPATH="mpif90" -D CMAKE_CXX_COMPILER:FILEPATH="mpicxx" -D CMAKE_C_COMPILER:FILEPATH="mpicc" \
-D BLAS_LIBRARY_DIRS="/home/issp/intel/composer_xe_2013/mkl/lib/intel64/"  -D TPL_BLAS_LIBRARIES:STRING="-mkl=parallel" \
-D LAPACK_LIBRARY_DIRS="/home/issp/intel/composer_xe_2013/mkl/lib/intel64/"  -D TPL_LAPACK_LIBRARIES:STRING="-mkl=parallel" \
-D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
-D TPL_ENABLE_Boost:BOOL=ON -D Boost_INCLUDE_DIRS:PATH=/opt/nano/alps/alps-20121208-r6630/include/boost/ -D Boost_LIBRARY_DIRS:PATH=/opt/nano/alps/alps-20121208-r6630/lib/ \
-D Trilinos_ENABLE_Anasazi:BOOL=ON \
-D Trilinos_ENABLE_Didasko:BOOL=ON \
-D Trilinos_ENABLE_EXAMPLES:BOOL=ON -D Trilinos_ENABLE_TESTS:BOOL=ON \
$HOME/build/trilinos-11.4.1-Source

make -j4
make test
make install
