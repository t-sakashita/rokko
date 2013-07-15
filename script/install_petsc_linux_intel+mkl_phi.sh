#!/bin/bash -x

PREFIX="$1"
#test -z "$PREFIX" && PREFIX="${HOME}/lib/petsc-3.4.2_installed/${PETSC_ARCH}"
test -z "$PREFIX" && PREFIX="${HOME}/rokko_lib"
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf petsc-3.4.2
wget -O - http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.2.tar.gz | tar zxf -

cd petsc-3.4.2
unset PETSC_DIR

./configure --prefix="$PREFIX" --with-cc="mpicc" --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS=-O3 --with-cxx="mpicxx" --with-fc="mpif90" --with-mpiexec="mpiexec" --with-blas-lapack-dir=/home/issp/intel/composer_xe_2013.4.183/mkl/bin/intel64 --with-c++-support=1 --with-debugging=0
make
make install
# You can't use -j4 option no longer in make, SLEPc automatically sets the number of parallelization.
