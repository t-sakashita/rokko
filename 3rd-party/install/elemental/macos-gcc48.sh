#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix
set_build_dir

. $PREFIX_OPT/env.sh

cd $BUILD_DIR
rm -rf elemental-$ELEMENTAL_VERSION
if [ -f $HOME/source/elemental-$ELEMENTAL_VERSION.tgz ]; then
  check tar zxf $HOME/source/elemental-$ELEMENTAL_VERSION.tgz
else
  check wget -O - http://elemental.googlecode.com/files/elemental-$ELEMENTAL_VERSION.tgz | tar zxf -
fi
cd $BUILD_DIR/elemental-$ELEMENTAL_VERSION
patch -p1 < $SCRIPT_DIR/elemental-$ELEMENTAL_VERSION.patch

BUILD_TYPES="Release Debug"

for build_type in $BUILD_TYPES; do
  cd $BUILD_DIR
  rm -rf elemental-$ELEMENTAL_VERSION-build-$build_type
  mkdir -p elemental-$ELEMENTAL_VERSION-build-$build_type && cd elemental-$ELEMENTAL_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE="Hybrid$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX_ROKKO/$build_type \
    -DCMAKE_C_COMPILER=gcc-mp-4.8 -DCMAKE_Fortran_COMPILER=gfortran-mp-4.8 \
    -DMATH_LIBS="-Wl,-framework -Wl,vecLib" \
    -DSHARED_LIBRARIES=ON \
    -DCMAKE_INSTALL_RPATH="$PREFIX_ROKKO/$build_type/lib" -DCMAKE_SKIP_BUILD_RPATH=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_MACOSX_RPATH=1 \
    $BUILD_DIR/elemental-$ELEMENTAL_VERSION
  check make VERBOSE=1 -j4
  $SUDO make install
done
