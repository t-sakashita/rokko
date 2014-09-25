#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix
set_build_dir

. $PREFIX_OPT/env.sh

cd $BUILD_DIR
rm -rf scalapack-$SCALAPACK_VERSION scalapack-$SCALAPACK_VERSION-build
if test -f $HOME/source/scalapack-$SCALAPACK_VERSION.tgz; then
  check tar zxf $HOME/source/scalapack-$SCALAPACK_VERSION.tgz
else
  check wget -O - http://www.netlib.org/scalapack/scalapack-$SCALAPACK_VERSION.tgz | tar zxf -
fi

BUILD_TYPES="Release Debug"

for build_type in $BUILD_TYPES; do
  cd $BUILD_DIR
  rm -rf scalapack-$SCALAPACK_VERSION-build-$build_type
  mkdir -p scalapack-$SCALAPACK_VERSION-build-$build_type
  cd scalapack-$SCALAPACK_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX_ROKKO/$build_type \
    -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
    -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
    -DLAPACK_LIBRARIES="-Wl,-framework -Wl,vecLib" \
    -DCMAKE_INSTALL_NAME_DIR=$$PREFIX_ROKKO/$build_type/lib \
    $BUILD_DIR/scalapack-$SCALAPACK_VERSION
  check make VERBOSE=1 -j4
  $SUDO make install
done
