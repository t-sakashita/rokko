#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION/$build_type
  cd $BUILD_DIR
  rm -rf Elemental-$ELEMENTAL_VERSION-build-$build_type
  mkdir -p Elemental-$ELEMENTAL_VERSION-build-$build_type && cd Elemental-$ELEMENTAL_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE="Hybrid$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DELEM_SHARED_LIBRARIES=ON \
    $BUILD_DIR/Elemental-$ELEMENTAL_VERSION
  check make VERBOSE=1
  $SUDO make install
done

DATE=$(date +%Y%m%d-%H%M%S)
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION/$build_type
  cat << EOF > $BUILD_DIR/elementalvars.sh
# elemental $(basename $0 .sh) $ELEMENTAL_VERSION $ELEMENTAL_RK_REVISION $DATE
export ELEMENTAL_ROOT=$PREFIX_ROKKO/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION
export LD_LIBRARY_PATH=$PREFIX/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/elementalvars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/elementalvars.sh
# elemental $(basename $0 .sh) $ELEMENTAL_VERSION $ELEMENTAL_RK_REVISION $DATE
export ELEMENTAL_ROOT=$PREFIX_ROKKO/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION
EOF
$SUDO cp -f $BUILD_DIR/elementalvars.sh $PREFIX_ROKKO/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION
