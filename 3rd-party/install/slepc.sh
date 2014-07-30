#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf slepc-3.5.0
if test -f $HOME/source/slepc-3.5.0.tar.gz; then
  tar zxf $HOME/source/slepc-3.5.0.tar.gz
else
  wget -O - http://www.grycap.upv.es/slepc/download/download.php?filename=slepc-3.5.0.tar.gz | tar zxf -
fi

cd $HOME/build/slepc-3.5.0
rm -rf $PREFIX/include/slepc*
unset SLEPC_DIR PETSC_DIR PETSC_ARCH
./configure --prefix=$PREFIX

make SLEPC_DIR=$PWD PETSC_DIR=$PREFIX
make SLEPC_DIR=$PWD PETSC_DIR=$PREFIX install
