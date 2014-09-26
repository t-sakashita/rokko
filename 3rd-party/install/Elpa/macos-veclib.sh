#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix
set_build_dir

. $PREFIX_OPT/env.sh

cd $BUILD_DIR
rm -rf elpa_lib-$ELPA_VERSION
check tar zxf $HOME/source/elpa_lib-$ELPA_VERSION.tar.gz
cd $BUILD_DIR/elpa_lib-$ELPA_VERSION
patch -p1 < $SCRIPT_DIR/elpa_lib-201305.patch

BUILD_TYPES="Release Debug"

for build_type in $BUILD_TYPES; do
  cd $BUILD_DIR
  rm -rf elpa_lib-$ELPA_VERSION-build-$build_type
  mkdir -p elpa_lib-$ELPA_VERSION-build-$build_type
  cd elpa_lib-$ELPA_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX_ROKKO/$build_type \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_Fortran_FLAGS_RELEASE="-O3 -ffree-line-length-none" \
    -DCMAKE_Fortran_FLAGS_DEBUG="-g -ffree-line-length-none" \
    -DSCALAPACK_LIB="-L$PREFIX_ROKKO/$build_type/lib -lscalapack -Wl,-framework -Wl,vecLib" \
    -DCMAKE_INSTALL_RPATH="$PREFIX_ROKKO/$build_type/lib" -DCMAKE_SKIP_BUILD_RPATH=OFF -DCMAKE_BUILD_WITH_INSTALL_RPATH=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_MACOSX_RPATH=1 \
    $BUILD_DIR/elpa_lib-$ELPA_VERSION
  check make VERBOSE=1 -j4
  $SUDO make install
done
