#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix
set_build_dir

. $PREFIX_OPT/env.sh
. $PREFIX_OPT/gcc-4.9.sh

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  cd $BUILD_DIR
  rm -rf Elemental-$ELEMENTAL_VERSION-build-$build_type
  mkdir -p Elemental-$ELEMENTAL_VERSION-build-$build_type && cd Elemental-$ELEMENTAL_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE="Hybrid$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX_ROKKO/$build_type \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran \
    -DELEM_SHARED_LIBRARIES=ON \
    $BUILD_DIR/Elemental-$ELEMENTAL_VERSION
  check make VERBOSE=1 -j4
  $SUDO make install
done
