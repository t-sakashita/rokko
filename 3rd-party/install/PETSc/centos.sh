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
    ./configure --prefix=$PREFIX_ROKKO/$build_type \
      --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS="-O3" \
      --with-scalapack-lib="/usr/lib64/openmpi/lib/libscalapack.so /usr/lib64/openmpi/lib/libmpiblacs.so /usr/lib64/openmpi/lib/libmpiblacsF77init.so" \
      --with-c++-support=1 --with-debugging=0
  else
    ./configure --prefix=$PREFIX_ROKKO/$build_type \
      --COPTFLAGS="-g -O0" --CXXOPTFLAGS="-g -O0" --FOPTFLAGS="-g -O0" \
      --with-scalapack-lib="/usr/lib64/openmpi/lib/libscalapack.so /usr/lib64/openmpi/lib/libmpiblacs.so /usr/lib64/openmpi/lib/libmpiblacsF77init.so" \
      --with-c++-support=1 --with-debugging=1
  fi
  check make
  $SUDO make install
done
