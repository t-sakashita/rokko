#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf scalapack-2.0.2 scalapack-2.0.2-build
if test -f $HOME/source/scalapack-2.0.2.tgz; then
  tar zxf $HOME/source/scalapack-2.0.2.tgz
else
  wget -O - http://www.netlib.org/scalapack/scalapack-2.0.2.tgz | tar zxf -
fi

cd $HOME/build
mkdir -p scalapack-2.0.2-build && cd scalapack-2.0.2-build
cmake -DCMAKE_C_COMPILER=openmpicc -DCMAKE_Fortran_COMPILER=openmpif90 \
    -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
    -DLAPACK_LIBRARIES="-L/opt/local/lib -llapack -lptcblas -lptf77blas -latlas -lgfortran" \
    -DCMAKE_INSTALL_NAME_DIR=$PREFIX/lib \
    -DCMAKE_INSTALL_PREFIX=$PREFIX $HOME/build/scalapack-2.0.2

make -j4
make install
