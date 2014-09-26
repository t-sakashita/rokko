#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix
set_build_dir

. $PREFIX_OPT/env.sh

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  cd $BUILD_DIR
  mkdir -p elpa_lib-$ELPA_VERSION-build-$build_type
  cd elpa_lib-$ELPA_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX_ROKKO/$build_type \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_Fortran_FLAGS_RELEASE="-O3 -ffree-line-length-none" \
    -DCMAKE_Fortran_FLAGS_DEBUG="-g -ffree-line-length-none" \
    -DSCALAPACK_LIB="-L/usr/lib64/openmpi/lib -lscalapack -lmpiblacs -lmpiblacsF77init -llapack -lblas" \
    $BUILD_DIR/elpa_lib-$ELPA_VERSION
  check make VERBOSE=1 -j4
  $SUDO make install
done
