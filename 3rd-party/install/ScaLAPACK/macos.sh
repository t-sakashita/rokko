#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/scalapack-$SCALAPACK_VERSION-$SCALAPACK_PATCH_VERSION/$build_type
  cd $BUILD_DIR
  rm -rf scalapack-$SCALAPACK_VERSION-build-$build_type
  mkdir -p scalapack-$SCALAPACK_VERSION-build-$build_type
  cd scalapack-$SCALAPACK_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
    -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
    -DLAPACK_LIBRARIES="-Wl,-framework -Wl,Accelerate" \
    -DCMAKE_INSTALL_RPATH="$PREFIX/lib" -DCMAKE_SKIP_BUILD_RPATH=OFF -DCMAKE_BUILD_WITH_INSTALL_RPATH=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_MACOSX_RPATH=1 \
    $BUILD_DIR/scalapack-$SCALAPACK_VERSION
  check make VERBOSE=1 -j4
  $SUDO make install
  cat << EOF > $BUILD_DIR/scalapackvars.sh
export SCALAPACK_ROOT=$PREFIX
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/scalapackvars.sh $PREFIX
done
