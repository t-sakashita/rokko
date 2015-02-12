#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  $SUDO rm -f $PREFIX_ROKKO/slepc
  $SUDO ln -s $PREFIX_ROKKO/slepc-$SLEPC_VERSION-$SLEPC_PATCH_VERSION $PREFIX_ROKKO/slepc
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/slepcvars.sh
  $SUDO ln -s $PREFIX_ROKKO/slepc/$build_type/slepcvars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/slepcvars.sh
done

$SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/slepcvars.sh
$SUDO ln -s $PREFIX_ROKKO/slepc/slepcvars.sh $PREFIX_ROKKO/rokkoenv.d/slepcvars.sh
