#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX="/opt/nano/rokko"
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build && cd $HOME/build
rm -rf elemental-0.78-p1 elemental-0.78-p1-build
wget -O - http://elemental.googlecode.com/files/elemental-0.78-p1.tgz | tar zxf -
mkdir -p elemental-0.78-p1-build && cd elemental-0.78-p1-build
cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 -DMATH_LIBS="-mkl=parallel;-lifcore" -DIFCORE_LIB="-lifcore" -DSHARED_LIBRARIES=ON -DCMAKE_INSTALL_PREFIX=$PREFIX $HOME/build/elemental-0.78-p1

make -j4 all
make install
