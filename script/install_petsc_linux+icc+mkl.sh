#!/bin/bash -x

PREFIX="$1"
if test -z "$PREFIX"; then
    if test -d "/opt/nano/rokko"; then
        PREFIX=/opt/nano/rokko
    elif test -d "/opt/rokko"; then
        PREFIX=/opt/nano/rokko
    fi
fi
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf petsc-3.3-p6
wget -O - http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.3-p6.tar.gz | tar zxf -

cd petsc-3.3-p6
./configure --prefix="$PREFIX" --with-cc="mpicc" --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS=-O3 --with-cxx="mpic++" --with-fc="mpif77" --with-mpiexec="mpiexec" --with-blas-lapack-lib="-mkl=parallel" --with-blacs=1 --with-scalapack=1 --with-c++-support=1 --with-debugging=0
make -j4
make install
