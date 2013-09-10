#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"
dir=`dirname $0`
SCRIPT_DIR=`cd $dir && pwd`

mkdir -p $HOME/build
cd $HOME/build
rm -rf elemental-0.80 elemental-0.80-build
if test -f $HOME/source/elemental-0.80.tgz; then
  tar zxf $HOME/source/elemental-0.80.tgz
else
  wget -O - http://elemental.googlecode.com/files/elemental-0.80.tgz | tar zxf -
fi

cd $HOME/build/elemental-0.80
patch -p1 < $SCRIPT_DIR/elemental-0.80.patch

mkdir elemental-0.80-build && cd elemental-0.80-build
cmake -DCMAKE_CXX_COMPILER=openmpicxx -DCMAKE_C_COMPILER=openmpicc \
    -DCMAKE_Fortran_COMPILER=openmpif90 -DSHARED_LIBRARIES=ON \
    -DMATH_LIBS="-Wl,-framework -Wl,vecLib" \
    -DCMAKE_INSTALL_NAME_DIR=$PREFIX/lib -DCMAKE_INSTALL_PREFIX=$PREFIX -DELEM_EXAMPLES=ON \
    -DELEM_TESTS=ON $HOME/build/elemental-0.80

make -j4 all
make install
