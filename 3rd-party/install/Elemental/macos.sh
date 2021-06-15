#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/elemental/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION/$build_type
  cd $BUILD_DIR
  rm -rf elemental-$ELEMENTAL_VERSION-build-$build_type
  mkdir -p elemental-$ELEMENTAL_VERSION-build-$build_type && cd elemental-$ELEMENTAL_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE="Hybrid$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran \
    -DELEM_SHARED_LIBRARIES=ON \
    -DCMAKE_INSTALL_RPATH="$PREFIX/lib" -DCMAKE_SKIP_BUILD_RPATH=OFF -DCMAKE_BUILD_WITH_INSTALL_RPATH=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_MACOSX_RPATH=1 \
    -DEL_IGNORE_OSX_GCC_ALIGNMENT_PROBLEM=ON \
    $BUILD_DIR/elemental-$ELEMENTAL_VERSION
  check make -j4
  $SUDO make install
done

DATE=$(date +%Y%m%d-%H%M%S)
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/elemental/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION/$build_type
  cat << EOF > $BUILD_DIR/elementalvars.sh
# elemental $(basename $0 .sh) $ELEMENTAL_VERSION $ELEMENTAL_RK_REVISION $DATE
export ELEMENTAL_ROOT=$PREFIX_ROKKO/elemental/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION
export LD_LIBRARY_PATH=$PREFIX/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/elementalvars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/elementalvars.sh
# elemental $(basename $0 .sh) $ELEMENTAL_VERSION $ELEMENTAL_RK_REVISION $DATE
export ELEMENTAL_ROOT=$PREFIX_ROKKO/elemental/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION
EOF
$SUDO cp -f $BUILD_DIR/elementalvars.sh $PREFIX_ROKKO/elemental/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION
