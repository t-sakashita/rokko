#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix
set_build_dir

. $PREFIX_OPT/env.sh

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"

unset PETSC_DIR
for build_type in $BUILD_TYPES; do
  cd $BUILD_DIR
  cp -rp petsc-$PETSC_VERSION petsc-$PETSC_VERSION-build-$build_type
  cd petsc-$PETSC_VERSION-build-$build_type
  if [ $build_type == "Release" ]; then
    if [ `which mpicxx > /dev/null 2>&1; echo $?` = 0 ]; then
      ./configure --prefix=$PREFIX_ROKKO/$build_type \
        --with-cxx=mpicxx --with-cc=mpicc --with-fc=mpif90 \
        --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS="-O3" \
        --with-mpiexec="mpiexec" \
        --with-blas-lapack-dir=$MKLROOT/bin/intel64 \
        --with-c++-support=1 --with-debugging=0
    else
      ./configure --prefix=$PREFIX_ROKKO/$build_type \
        --with-cxx=icpc --with-cc=icc --with-fc=ifort \
        --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS="-O3" \
        --with-mpi-lib="-lmpi" \
        --with-blas-lapack-dir=$MKLROOT/bin/intel64 \
        --with-c++-support=1 --with-debugging=0
    fi
  else
    if [ `which mpicxx > /dev/null 2>&1; echo $?` = 0 ]; then
      ./configure --prefix=$PREFIX_ROKKO/$build_type \
        --with-cxx=mpicxx --with-cc=mpicc --with-fc=mpif90 \
        --COPTFLAGS="-O0 -g" --CXXOPTFLAGS="-O0 -g" --FOPTFLAGS="-O0 -g" \
        --with-mpiexec="mpiexec" \
        --with-blas-lapack-dir=$MKLROOT/bin/intel64 \
        --with-c++-support=1 --with-debugging=1
    else
      ./configure --prefix=$PREFIX_ROKKO/$build_type \
        --with-cxx=icpc --with-cc=icc --with-fc=ifort \
        --COPTFLAGS="-O0 -g" --CXXOPTFLAGS="-O0 -g" --FOPTFLAGS="-O0 -g" \
        --with-mpi-lib="-lmpi" \
        --with-blas-lapack-dir=$MKLROOT/bin/intel64 \
        --with-c++-support=1 --with-debugging=1
    fi
  fi
  check make
  $SUDO LD_LIBRARY_PATH=$LD_LIBRARY_PATH make install
done
