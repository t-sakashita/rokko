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
patch -p1 < $SCRIPT_DIR/elemental-0.80-fx10.patch

cd $HOME/build
mkdir elemental-0.80-build && cd elemental-0.80-build
cmake -DCMAKE_C_COMPILER=mpifccpx -DCMAKE_C_FLAGS="-Kfast -Xg -mt" \
    -DCMAKE_CXX_COMPILER=mpiFCCpx -DCMAKE_CXX_FLAGS="-Kfast -Xg -mt" \
    -DCMAKE_Fortran_COMPILER=mpifrtpx -DCMAKE_Fortran_FLAGS="-Kfast -mt" \
    -DMATH_LIBS="-SSL2 --linkfortran" \
    -DSHARED_LIBRARIES=ON \
    -DCMAKE_INSTALL_PREFIX=$PREFIX $HOME/build/elemental-0.80

make -j4 all VERBOSE=1
make install
