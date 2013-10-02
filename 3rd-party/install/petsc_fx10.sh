#!/bin/bash

# Reference: http://www.openpetascale.org/source/PETSc%20Installation%20Guide%20(June%202013).pdf

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"
dir=`dirname $0`
SCRIPT_DIR=`cd $dir && pwd`

mkdir -p $HOME/build
cd $HOME/build
rm -rf petsc-3.4.2
if test -f $HOME/source/petsc-3.4.2.tar.gz; then
  tar zxf $HOME/source/petsc-3.4.2.tar.gz
else
  wget -O - http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.2.tar.gz | tar zxf -
fi

cd $HOME/build/petsc-3.4.2
patch -p0 < $SCRIPT_DIR/petsc-3.4.2-fx10.patch

cd $HOME/build/petsc-3.4.2
./configure --prefix=$PREFIX \
  --with-cc="mpifccpx" --CFLAGS="-mt -Xg" --COPTFLAGS="-Kfast" \
  --with-cxx="mpiFCCpx" --CXXFLAGS="-mt -Xg" --CXXOPTFLAGS="-Kfast" \
  --with-fc="mpifrtpx" --FFLAGS="-Kthreadsafe" --FOPTFLAGS="-Kfast" \
  --with-blas-lapack-lib="-SSL2" \
  --with-x=0 --with-c++-support=1 --with-info=1 --with-debugging=0 --known-mpi-shared-libraries=0 --with-valgrind=0 \
  --with-batch=1 \
  --known-level1-dcache-size=32768 \
  --known-level1-dcache-linesize=32 \
  --known-level1-dcache-assoc=0 \
  --known-memcmp-ok=1 \
  --known-sizeof-char=1 \
  --known-sizeof-void-p=8 \
  --known-sizeof-short=2 \
  --known-sizeof-int=4 \
  --known-sizeof-long=8 \
  --known-sizeof-long-long=8 \
  --known-sizeof-float=4 \
  --known-sizeof-double=8 \
  --known-sizeof-size_t=8 \
  --known-bits-per-byte=8 \
  --known-sizeof-MPI_Comm=8 \
  --known-sizeof-MPI_Fint=4 \
  --known-mpi-long-double=1

make # -j option can not be specified
make install
