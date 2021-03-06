#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_RK_REVISION/$build_type
  cd $BUILD_DIR
  mkdir -p trilinos-$TRILINOS_VERSION-Source-build-$build_type
  cd trilinos-$TRILINOS_VERSION-Source-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE="$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
    -DTPL_ENABLE_MPI=ON -DMPI_USE_COMPILER_WRAPPERS=ON \
    -DTrilinos_ENABLE_Anasazi=ON -DTrilinos_ENABLE_Didasko=ON -DTrilinos_ENABLE_EXAMPLES=ON -DTrilinos_ENABLE_TESTS=ON \
    $BUILD_DIR/trilinos-$TRILINOS_VERSION-Source
  check make VERBOSE=1
  $SUDO make install
done

DATE=$(date +%Y%m%d-%H%M%S)
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_RK_REVISION/$build_type
  cat << EOF > $BUILD_DIR/trilinosvars.sh
# trilinos $(basename $0 .sh) $TRILINOS_VERSION $TRILINOS_RK_REVISION $DATE
export TRILINOS_ROOT=$PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_RK_REVISION
export LD_LIBRARY_PATH=$PREFIX/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/trilinosvars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/trilinosvars.sh
# trilinos $(basename $0 .sh) $TRILINOS_VERSION $TRILINOS_RK_REVISION $DATE
export TRILINOS_ROOT=$PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_RK_REVISION
EOF
$SUDO cp -f $BUILD_DIR/trilinosvars.sh $PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_RK_REVISION
