#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  $SUDO rm -f $PREFIX_ROKKO/scalapack
  $SUDO ln -s $PREFIX_ROKKO/scalapack-$SCALAPACK_VERSION-$SCALAPACK_PATCH_VERSION $PREFIX_ROKKO/scalapack
  $SUDO rm -f $PREFIX_ROKKO/env-$build_type.d/scalapack.sh
  $SUDO ln -s $PREFIX_ROKKO/scalapack/$build_type/scalapackvars.sh $PREFIX_ROKKO/env-$build_type.d/scalapackvars.sh
done
