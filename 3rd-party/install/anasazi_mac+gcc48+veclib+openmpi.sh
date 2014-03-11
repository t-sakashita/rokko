#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf trilinos-11.4.2-Source trilinos-11.4.2-Source-build
tar jxf $HOME/source/trilinos-11.4.2-Source.tar.bz2

cd $HOME/build
mkdir -p trilinos-11.4.2-Source-build && cd trilinos-11.4.2-Source-build
cmake \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
    -DTPL_ENABLE_MPI=ON -DMPI_USE_COMPILER_WRAPPERS=ON \
    -DCMAKE_Fortran_COMPILER="mpif90-openmpi-gcc48" -DCMAKE_CXX_COMPILER="mpicxx-openmpi-gcc48" -DCMAKE_C_COMPILER="mpicc-openmpi-gcc48" \
    -DTPL_BLAS_LIBRARIES="-framework veclib" \
    -DTPL_LAPACK_LIBRARIES="-framework veclib" \
    -DTrilinos_ENABLE_Anasazi=ON \
    -DTrilinos_ENABLE_Didasko=ON \
    -DTrilinos_ENABLE_EXAMPLES=ON -DTrilinos_ENABLE_TESTS=ON \
    $HOME/build/trilinos-11.4.2-Source

make -j4 VERBOSE=1
make test
make install
