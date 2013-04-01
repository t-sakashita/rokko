#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf slepc-3.3-p3
wget -O - http://www.grycap.upv.es/slepc/download/download.php?filename=slepc-3.3-p3.tar.gz | tar zxf -

cd slepc-3.3-p3
export PETSC_DIR=$PREFIX
./configure --prefix=$PREFIX

export SLEPC_DIR=$PWD
export PETSC_ARCH=arch-installed-petsc
make -j4
make install
