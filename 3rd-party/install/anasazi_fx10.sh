#!/bin/bash

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

TMP_DIR=`dirname $0`
SCRIPT_DIR=`cd $TMP_DIR && pwd`

mkdir -p $HOME/build
cd $HOME/build
rm -rf trilinos-11.8.1-Source trilinos-11.8.1-Source-build
tar jxf $HOME/source/trilinos-11.8.1-Source.tar.bz2

cd $HOME/build/trilinos-11.8.1-Source
patch -p1 < $SCRIPT_DIR/trilinos-11.8.1-Source.patch

cd $HOME/build
mkdir -p trilinos-11.8.1-Source-build && cd trilinos-11.8.1-Source-build
cmake \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
  -DTPL_ENABLE_MPI=ON \
  -DCMAKE_Fortran_COMPILER=mpifrtpx -DCMAKE_CXX_COMPILER=mpiFCCpx -DCMAKE_C_COMPILER=mpifccpx \
  -DCMAKE_CXX_FLAGS="-Kfast -KPIC -Kopenmp -Xg -mt" -DCMAKE_C_FLAGS="-Kfast -KPIC -Kopenmp -Xg -mt" \
  -DCMAKE_Fortran_FLAGS="-Kfast -KPIC -Cpp"\
  -DTPL_BLAS_LIBRARIES="-SSL2BLAMP --linkfortran" \
  -DTPL_LAPACK_LIBRARIES="-SSL2BLAMP --linkfortran" \
  -DCMAKE_EXE_LINKER_FLAGS="--linkfortran" \
  -DTrilinos_EXTRA_LINK_FLAGS="--linkfortran" \
  -DTrilinos_ENABLE_Anasazi=ON \
  -DTrilinos_ENABLE_Didasko=ON \
  -DTrilinos_ENABLE_EXAMPLES=ON -D Trilinos_ENABLE_TESTS=ON \
  -DTrilinos_SKIP_FORTRANCINTERFACE_VERIFY_TEST=ON \
$HOME/build/trilinos-11.8.1-Source

make -j4 VERBOSE=1
make test
make install
