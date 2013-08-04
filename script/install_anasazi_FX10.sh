#!/bin/bash

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

TMP_DIR=`dirname $0`
SCRIPT_DIR=`cd $TMP_DIR && pwd`

mkdir -p $WORK/build
cd $WORK/build
rm -rf trilinos-11.4.1-Source
tar jxf $WORK/source/trilinos-11.4.1-Source.tar.bz2

cd $WORK/build/trilinos-11.4.1-Source
patch -p1 < $SCRIPT_DIR/trilinos-11.4.1-Source.patch

cd $WORK/build
rm -rf trilinos-11.4.1-build
mkdir -p trilinos-11.4.1-build
cd trilinos-11.4.1-build

cmake \
-D TPL_ENABLE_MPI:BOOL=ON  \
-D CMAKE_Fortran_COMPILER:FILEPATH="mpifrtpx" -D CMAKE_CXX_COMPILER:FILEPATH="mpiFCCpx" -D CMAKE_C_COMPILER:FILEPATH="mpifccpx" \
-D CMAKE_CXX_FLAGS="-Kopenmp -Xg -mt" \
-D CMAKE_C_FLAGS="-Kopenmp -Xg -mt" \
-D TPL_BLAS_LIBRARIES:STRING="-SSL2BLAMP --linkfortran" \
-D TPL_LAPACK_LIBRARIES:STRING="-SSL2BLAMP --linkfortran" \
-D CMAKE_INSTALL_PREFIX:PATH="../install_build/" \
-D CMAKE_EXE_LINKER_FLAGS="--linkfortran" \
-D Trilinos_EXTRA_LINK_FLAGS:STRING="--linkfortran" \
-D Trilinos_ENABLE_Anasazi:BOOL=ON \
-D Trilinos_ENABLE_EXAMPLES:BOOL=ON -D Trilinos_ENABLE_TESTS:BOOL=ON \
-D Trilinos_SKIP_FORTRANCINTERFACE_VERIFY_TEST=ON \
$WORK/build/trilinos-11.4.1-Source

make -j8 2>&1 | tee make.log
make install
