#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh
cd $BUILD_DIR/petsc-$PETSC_VERSION
patch -p1 < $SCRIPT_DIR/petsc-3.5.2-fx10.patch

unset PETSC_DIR
BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX_BACKEND=$PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION/Linux-s64fx/$build_type
  cd $BUILD_DIR
  cp -rp petsc-$PETSC_VERSION petsc-$PETSC_VERSION-build-Linux-s64fx-$build_type
  cd petsc-$PETSC_VERSION-build-Linux-s64fx-$build_type
  if [ $build_type == "Release" ]; then
    ./configure --prefix=$PREFIX_BACKEND \
      --with-cc="mpifccpx" --CFLAGS="-mt -Xg" --COPTFLAGS="-Kfast" \
      --with-cxx="mpiFCCpx" --CXXFLAGS="-mt -Xg" --CXXOPTFLAGS="-Kfast" \
      --with-fc="mpifrtpx" --FFLAGS="-Kthreadsafe" --FOPTFLAGS="-Kfast" \
      --LDFLAGS="-lmpi_f77 -lmpi_f90" \
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
  else
    ./configure --prefix=$PREFIX_BACKEND \
      --with-cc="mpifccpx" --CFLAGS="-mt -Xg" --COPTFLAGS="-Kfast" \
      --with-cxx="mpiFCCpx" --CXXFLAGS="-mt -Xg" --CXXOPTFLAGS="-Kfast" \
      --with-fc="mpifrtpx" --FFLAGS="-Kthreadsafe" --FOPTFLAGS="-Kfast" \
      --LDFLAGS="-lmpi_f77 -lmpi_f90" \
      --with-blas-lapack-lib="-SSL2" \
      --with-x=0 --with-c++-support=1 --with-info=1 --with-debugging=1 --known-mpi-shared-libraries=0 --with-valgrind=0 \
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
  fi
  check make
  $SUDO make install
done
