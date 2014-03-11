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

cd $HOME/build
mkdir elemental-0.80-build && cd elemental-0.80-build
cmake -DCMAKE_CXX_COMPILER=mpicxx-openmpi-gcc48 -DCMAKE_C_COMPILER=mpicc-openmpi-gcc48 \
    -DCMAKE_Fortran_COMPILER=mpif90-openmpi-gcc48 \
    -DMATH_LIBS="-Wl,-framework -Wl,vecLib" \
    -DCMAKE_INSTALL_NAME_DIR=$PREFIX/lib \
    -DELEM_EXAMPLES=ON -DELEM_TESTS=ON \
    -DSHARED_LIBRARIES=ON \
    -DCMAKE_INSTALL_PREFIX=$PREFIX $HOME/build/elemental-0.80

make -j4 all VERBOSE=1
make install
