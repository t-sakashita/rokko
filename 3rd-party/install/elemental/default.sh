#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix
set_build_dir

. $PREFIX_OPT/env.sh

cd $BUILD_DIR
rm -rf Elemental-$ELEMENTAL_VERSION
if [ -f $HOME/source/Elemental-$ELEMENTAL_VERSION.tgz ]; then
  check tar zxf $HOME/source/Elemental-$ELEMENTAL_VERSION.tgz
else
  check wget -O - http://elemental.googlecode.com/files/Elemental-$ELEMENTAL_VERSION.tgz | tar zxf -
fi
cd $BUILD_DIR/Elemental-$ELEMENTAL_VERSION
patch -p1 < $SCRIPT_DIR/elemental-$ELEMENTAL_VERSION.patch

BUILD_TYPES="Release Debug"

for build_type in $BUILD_TYPES; do
  cd $BUILD_DIR
  rm -rf Elemental-$ELEMENTAL_VERSION-build-$build_type
  mkdir -p Elemental-$ELEMENTAL_VERSION-build-$build_type && cd Elemental-$ELEMENTAL_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE="Hybrid$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX_ROKKO/$build_type \
    -DSHARED_LIBRARIES=ON \
    $BUILD_DIR/Elemental-$ELEMENTAL_VERSION
  check make VERBOSE=1 -j4
  $SUDO make install
done

