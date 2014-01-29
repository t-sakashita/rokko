#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"
dir=`dirname $0`
SCRIPT_DIR=`cd $dir && pwd`

mkdir -p $HOME/build
cd $HOME/build
rm -rf elpa_lib-201305 elpa_lib-201305-build
tar zxf $HOME/source/elpa_lib-201305.tar.gz

cd $HOME/build/elpa_lib-201305
patch -p1 < $SCRIPT_DIR/elpa_lib-201305.patch

cd $HOME/build
mkdir -p elpa_lib-201305-build && cd elpa_lib-201305-build
cmake -DCMAKE_C_COMPILER=openmpicc \
    -DCMAKE_Fortran_COMPILER=openmpif90 -DCMAKE_Fortran_FLAGS="-O3 -ffree-line-length-none" \
    -DSCALAPACK_LIB="-L$PREFIX/lib -lscalapack -Wl,-framework -Wl,vecLib" \
    -DCMAKE_INSTALL_PREFIX=$PREFIX $HOME/build/elpa_lib-201305

make -j4 VERBOSE=1
make install
