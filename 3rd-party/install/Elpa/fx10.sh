#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_PATCH_VERSION/$build_type
  cd $BUILD_DIR
  mkdir -p elpa_lib-$ELPA_VERSION-build-$build_type
  cd elpa_lib-$ELPA_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_C_COMPILER=mpifccpx -DCMAKE_Fortran_COMPILER=mpifrtpx \
    -DCMAKE_Fortran_FLAGS_RELEASE="-Kfast" \
    -DCMAKE_Fortran_FLAGS_DEBUG="-g -O0" \
    -DSCALAPACK_LIB="-SCALAPACK -SSL2" \
    $BUILD_DIR/elpa_lib-$ELPA_VERSION
  check make VERBOSE=1 -j4
  $SUDO make install
  cat << EOF > $BUILD_DIR/elpavars.sh
export ELPA_ROOT=$PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_PATCH_VERSION
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/elpavars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/elpavars.sh
export ELPA_ROOT=$PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_PATCH_VERSION
EOF
$SUDO cp -f $BUILD_DIR/elpavars.sh $PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_PATCH_VERSION
