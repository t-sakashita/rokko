#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"
dir=`dirname $0`
SCRIPT_DIR=`cd $dir && pwd`

mkdir -p $HOME/build
cd $HOME/build
rm -rf EigenExa-1.0-RC5
tar zxf $HOME/source/EigenExa-1.0-RC5.tgz

cd $HOME/build/EigenExa-1.0-RC5
patch -p1 < $SCRIPT_DIR/EigenExa-1.0-RC5.patch

cd $HOME/build
rm -rf EigenExa-1.0-RC5-build && mkdir -p EigenExa-1.0-RC5-build && cd EigenExa-1.0-RC5-build
cmake -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_Fortran_FLAGS="-g" -DSCALAPACK_LIB="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" $HOME/build/EigenExa-1.0-RC5

make -j4
make install
