#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX_BACKEND=$PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_RK_REVISION/Linux-s64fx/$build_type
  cd $BUILD_DIR
  mkdir -p elpa_lib-$ELPA_VERSION-build-Linux-s64fx-$build_type
  cd elpa_lib-$ELPA_VERSION-build-Linux-s64fx-$build_type
  check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX_BACKEND \
    -DCMAKE_C_COMPILER=mpifccpx -DCMAKE_Fortran_COMPILER=mpifrtpx \
    -DCMAKE_Fortran_FLAGS_RELEASE="-Kfast" \
    -DCMAKE_Fortran_FLAGS_DEBUG="-g -O0" \
    -DSCALAPACK_LIB="-SCALAPACK -SSL2" \
    $BUILD_DIR/elpa_lib-$ELPA_VERSION
  check make VERBOSE=1 -j4
  $SUDO make install
done

DATE=$(date +%Y%m%d-%H%M%S)
for build_type in $BUILD_TYPES; do
  PREFIX_BACKEND=$PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_RK_REVISION/Linux-s64fx/$build_type
  cat << EOF > $BUILD_DIR/elpavars.sh
# elpa $(basename $0 .sh) $ELPA_VERSION $ELPA_RK_REVISION $DATE
export ELPA_ROOT=$PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_RK_REVISION/Linux-s64fx
export LD_LIBRARY_PATH=$PREFIX_BACKEND/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/elpavars.sh $PREFIX_BACKEND
done

cat << EOF > $BUILD_DIR/elpavars.sh
# elpa $(basename $0 .sh) $ELPA_VERSION $ELPA_RK_REVISION $DATE
export ELPA_ROOT=$PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_RK_REVISION/Linux-s64fx
EOF
$SUDO cp -f $BUILD_DIR/elpavars.sh $PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_RK_REVISION/Linux-s64fx
