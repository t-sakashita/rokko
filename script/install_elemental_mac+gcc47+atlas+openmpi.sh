#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf elemental-0.78-p1 elemental-0.78-p1-build
wget -O - http://elemental.googlecode.com/files/elemental-0.78-p1.tgz | tar zxf -
mkdir elemental-0.78-p1-build && cd elemental-0.78-p1-build
cmake -DCMAKE_CXX_COMPILER=openmpicxx -DCMAKE_C_COMPILER=openmpicc -DCMAKE_Fortran_COMPILER=openmpif90 -DMATH_LIBS="-L/opt/local/lib -llapack -lptcblas -lptf77blas -latlas -lgfortran" -DCMAKE_INSTALL_PREFIX=$PREFIX -DELEM_EXAMPLES=ON -DELEM_TESTS=ON $HOME/build/elemental-0.78-p1

make -j4 all
make install
