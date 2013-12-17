#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"
DEBUG="$2"
test -z "$DEBUG" && DEBUG=0
echo "DEBUG = $DEBUG"

mkdir -p $HOME/build
cd $HOME/build
rm -rf petsc-3.4.3
if test -f $HOME/source/petsc-3.4.3.tar.gz; then
  tar zxf $HOME/source/petsc-3.4.3.tar.gz
else
  wget -O - http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.3.tar.gz | tar zxf -
fi

cd $HOME/build/petsc-3.4.3
if [ "$DEBUG" = "0" ]; then
    OPTFLAGS="-O3 -xSSE3"
else
    OPTFLAGS="-O0 -g"
fi
unset PETSC_DIR
if [ `which mpicxx > /dev/null 2>&1; echo $?` = 0 ]; then
  ./configure --prefix="$PREFIX" \
      --with-cxx=mpicxx --with-cc=mpicc --with-fc=mpif90 \
      --CXXOPTFLAGS="$OPTFLAGS" --COPTFLAGS="$OPTFLAGS" --FOPTFLAGS="$OPTFLAGS" \
      --with-mpiexec="mpiexec" \
      --with-blas-lapack-dir=$MKLROOT/bin/intel64 --with-c++-support=1 --with-debugging="$DEBUG"
else
  ./configure --prefix="$PREFIX" \
      --with-cxx=icpc --with-cc=icc --with-fc=ifort \
      --CXXOPTFLAGS="$OPTFLAGS" --COPTFLAGS="$OPTFLAGS" --FOPTFLAGS="$OPTFLAGS" \
      --with-mpi-lib="-lmpi" \
      --with-blas-lapack-dir=$MKLROOT/bin/intel64 --with-c++-support=1 --with-debugging="$DEBUG"
fi
make # -j option can not be specified
make install
