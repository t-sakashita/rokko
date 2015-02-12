#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  $SUDO rm -f $PREFIX_ROKKO/petsc
  $SUDO ln -s $PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION $PREFIX_ROKKO/petsc
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/petscvars.sh
  $SUDO ln -s $PREFIX_ROKKO/petsc/$build_type/petscvars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/petscvars.sh
done

$SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/petscvars.sh
$SUDO ln -s $PREFIX_ROKKO/petsc/petscvars.sh $PREFIX_ROKKO/rokkoenv.d/petscvars.sh
