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
    -DCMAKE_C_COMPILER=openmpicc -DCMAKE_Fortran_COMPILER=openmpif90 \
    -DSCALAPACK_LIB="-L$PREFIX/lib -lscalapack -Wl,-framework -Wl,vecLib" \
    $HOME/build/EigenExa-1.3

make -j4 VERBOSE=1
make install
