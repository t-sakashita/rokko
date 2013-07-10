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
cmake -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_CXX_COMPILER=openmpicxx -DCMAKE_C_COMPILER=openmpicc -DCMAKE_Fortran_COMPILER=openmpif90 -DSCALAPACK_LIB="-L$PREFIX/lib -lscalapack -Wl,-framework -Wl,vecLib" $HOME/build/EigenExa-1.0-RC5

make -j4
make install
