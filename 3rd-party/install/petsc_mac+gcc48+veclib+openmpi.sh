#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf petsc-3.4.3
if test -f $HOME/source/petsc-3.4.3.tar.gz; then
  tar zxf $HOME/source/petsc-3.4.3.tar.gz
else
  wget -O - http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.3.tar.gz | tar zxf -
fi

cd $HOME/build/petsc-3.4.3
unset PETSC_DIR
./configure --prefix=$PREFIX --with-cc=mpicc-openmpi-gcc48 --COPTFLAGS="-O3" --with-cxx=mpicxx-openmpi-gcc48 --CXXOPTFLAGS="-O3" --with-fc=mpif77-openmpi-gcc48 --FOPTFLAGS="-O3" --with-mpiexec=mpiexec-openmpi-gcc48 \
    --with-blas-lapack-lib="-framework vecLib" \
    --with-blacs-dir=$PREFIX --with-scalapack-dir=$PREFIX --with-c++-support=1 --with-debugging=0

make # -j option can not be specified
make install
