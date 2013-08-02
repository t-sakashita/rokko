#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX="$HOME/rokko_lib"
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf slepc-3.4.1
wget -O - http://www.grycap.upv.es/slepc/download/download.php?filename=slepc-3.4.1.tar.gz | tar zxf -

cd slepc-3.4.1
export SLEPC_DIR=$PWD
export PETSC_DIR=$PREFIX  #$HOME/build/petsc-3.4.2  #  /arch-linux2-c-opt
unset PETSC_ARCH

./configure --prefix=$PREFIX
make PETSC_ARCH=arch-installed-petsc
make PETSC_ARCH=arch-installed-petsc install
# You can't use -j4 option no longer in make, SLEPc automatically sets the number of parallelization.

