#!/bin/bash -x

#export PETSC_ARCH=arch-installed-petsc

PREFIX="$1"
#test -z "$PREFIX" && PREFIX="${HOME}/lib/petsc-3.3-p6_installed/${PETSC_ARCH}"
test -z "$PREFIX" && PREFIX="${HOME}/lib/petsc-3.3-p6_installed"
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf petsc-3.3-p6
wget -O - http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.3-p6.tar.gz | tar zxf -

cd petsc-3.3-p6
#export PETSC_DIR=
unset PETSC_DIR

./configure --prefix="$PREFIX" --with-cc="mpicc" --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS=-O3 --with-cxx="mpicxx" --with-fc="mpif90" --with-mpiexec="mpiexec" --with-blas-lapack-lib="-mkl=sequential" --with-c++-support=1 --with-debugging=0
make -j4
make install
