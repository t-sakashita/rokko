#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"
dir=`dirname $0`
SCRIPT_DIR=`cd $dir && pwd`

mkdir -p $HOME/build
cd $HOME/build
rm -rf EigenExa-1.3
if test -f $HOME/source/EigenExa-1.3.tgz; then
  tar zxf $HOME/source/EigenExa-1.3.tgz
else
  wget -O - http://www.aics.riken.jp/labs/lpnctrt/EigenExa-1.3.tgz | tar zxf -
fi
cd $HOME/build/EigenExa-1.3
patch -p1 < $SCRIPT_DIR/EigenExa-1.3.patch

cd $HOME/build
rm -rf EigenExa-1.3-build && mkdir -p EigenExa-1.3-build && cd EigenExa-1.3-build
cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
    -DCMAKE_CXX_FLAGS="-O3 -xSSE3" -DCMAKE_C_FLAGS="-O3 -xSSE3" -DCMAKE_Fortran_FLAGS="-O3 -xSSE3" \
    -DSCALAPACK_LIB="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" \
    $HOME/build/EigenExa-1.3

make -j4 VERBOSE=1
make install
