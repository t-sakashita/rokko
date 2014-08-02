#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf petsc-3.5.1
if test -f $HOME/source/petsc-3.5.1.tar.gz; then
  tar zxf $HOME/source/petsc-3.5.1.tar.gz
else
  wget -O - http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.5.1.tar.gz | tar zxf -
fi

cd $HOME/build/petsc-3.5.1
#unset PETSC_DIR
PETSC_DIR=$PWD

./configure --prefix=$PREFIX --with-cc=mpicc --COPTFLAGS="-O3" --with-cxx=mpicxx --CXXOPTFLAGS="-O3" --with-fc=mpif77 --FOPTFLAGS="-O3" --with-mpiexec=mpirun \
    --with-blas-lapack-lib="-framework vecLib" \
    --with-blacs-dir=$PREFIX --with-c++-support=1 --with-debugging=1 --with-valgrind-include=/usr/local/include/

make # -j option can not be specified
make install
