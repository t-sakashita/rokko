#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf scalapack-2.0.2
if test -f $HOME/source/scalapack-2.0.2.tgz; then
  tar zxf $HOME/source/scalapack-2.0.2.tgz
else
  wget -O - http://www.netlib.org/scalapack/scalapack-2.0.2 | tar zxf -
fi

cd scalapack-2.0.2
cmake -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_C_COMPILER=openmpicc -DCMAKE_Fortran_COMPILER=openmpif90 -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF $HOME/build/scalapack-2.0.2
make -j2
make install

# ln -s libscalapack.a $PREFIX/lib/libblacs.a
