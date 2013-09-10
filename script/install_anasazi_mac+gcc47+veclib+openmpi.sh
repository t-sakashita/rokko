#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf trilinos-11.4.1-Source
tar jxf $HOME/source/trilinos-11.4.1-Source.tar.bz2

cd $HOME/build
mkdir -p trilinos-11.4.1-Source-build && cd trilinos-11.4.1-Source-build
cmake \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DTPL_ENABLE_MPI=ON -DMPI_USE_COMPILER_WRAPPERS=ON \
    -DCMAKE_Fortran_COMPILER="openmpif90" -DCMAKE_CXX_COMPILER="openmpicxx" -DCMAKE_C_COMPILER="openmpicc" \
    -DTPL_BLAS_LIBRARIES="-framework veclib" \
    -DTPL_LAPACK_LIBRARIES="-framework veclib" \
    -DTrilinos_ENABLE_Anasazi=ON \
    -DTrilinos_ENABLE_Didasko=ON \
    -DTrilinos_ENABLE_EXAMPLES=ON -DTrilinos_ENABLE_TESTS=ON \
    $HOME/build/trilinos-11.4.1-Source

# -DTPL_ENABLE_Boost=ON -DBoost_INCLUDE_DIRS=/opt/boost/boost_1_53_0

make -j4
make test
make install
