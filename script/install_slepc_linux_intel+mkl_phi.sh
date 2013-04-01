#!/bin/bash -x

#export PETSC_ARCH=arch-installed-petsc

PREFIX="$1"
#test -z "$PREFIX" && PREFIX="${HOME}/lib/slepc-3.3-p3_installed/${PETSC_ARCH}"
test -z "$PREFIX" && PREFIX="${HOME}/lib/petsc-3.3-p6_installed"
echo "PREFIX = $PREFIX"


mkdir -p $HOME/build
cd $HOME/build
rm -rf slepc-3.3-p3
wget -O - http://www.grycap.upv.es/slepc/download/download.php?filename=slepc-3.3-p3.tar.gz | tar zxf -

cd slepc-3.3-p3
#export PETSC_DIR="${HOME}/lib/petsc-3.3-p6_installed/"   #${PETSC_ARCH}"  
export PETSC_DIR=$PREFIX

./configure --prefix=$PREFIX

export SLEPC_DIR=$PWD
export PETSC_ARCH=arch-installed-petsc
make -j4
make install
