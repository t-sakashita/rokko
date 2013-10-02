#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf petsc-3.4.2
if test -f $HOME/source/petsc-3.4.2.tar.gz; then
  tar zxf $HOME/source/petsc-3.4.2.tar.gz
else
  wget -O - http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.2.tar.gz | tar zxf -
fi

cd $HOME/build/petsc-3.4.2
./configure --prefix=$PREFIX --with-cc=openmpicc --COPTFLAGS="-O3" --with-cxx=openmpicxx --CXXOPTFLAGS="-O3" --with-fc=openmpif77 --FOPTFLAGS="-O3" --with-mpiexec=openmpiexec \
    --with-blas-lapack-lib="/opt/local/lib/liblapack.a -lptcblas -lptf77blas -latlas" \
    --with-blacs-dir=$PREFIX --with-scalapack-dir=$PREFIX --with-c++-support=1 --with-debugging=0

make # -j option can not be specified
make install
