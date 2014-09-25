#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"
dir=`dirname $0`
SCRIPT_DIR=`cd $dir && pwd`

mkdir -p $HOME/build
cd $HOME/build
rm -rf elemental-0.84-p1 elemental-0.84-p1-build
if test -f $HOME/source/elemental-0.84-p1.tgz; then
  tar zxf $HOME/source/elemental-0.84-p1.tgz
else
  wget -O - http://libelemental.org/pub/releases/Elemental-0.84-p1.tgz | tar zxf -
fi

cd $HOME/build/elemental-0.84-p1
patch -p1 < $SCRIPT_DIR/elemental-0.84-p1.patch

cd $HOME/build
mkdir elemental-0.84-p1-build && cd elemental-0.84-p1-build
cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_Fortran_COMPILER=mpif90 \
    -DMATH_LIBS="-Wl,-framework -Wl,vecLib" \
    -DCMAKE_INSTALL_NAME_DIR=$PREFIX/lib \
    -DELEM_EXAMPLES=ON -DELEM_TESTS=ON \
    -DSHARED_LIBRARIES=ON \
    -DCMAKE_INSTALL_PREFIX=$PREFIX $HOME/build/elemental-0.84-p1

make -j4 all VERBOSE=1
make install
