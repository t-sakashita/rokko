#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix
set_build_dir

. $PREFIX_OPT/env.sh

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  cd $BUILD_DIR
  mkdir -p EigenExa-$EIGENEXA_VERSION-build-$build_type
  cd EigenExa-$EIGENEXA_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX_ROKKO/$build_type \
    -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
    -DSCALAPACK_LIB="-L$PREFIX_ROKKO/$build_type/lib -L/opt/local/lib -lscalapack -Wl,-framework -Wl,vecLib" \
    $BUILD_DIR/EigenExa-$EIGENEXA_VERSION
  check make VERBOSE=1 -j4
  $SUDO make install
done
