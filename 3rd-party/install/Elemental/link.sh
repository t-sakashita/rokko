#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  $SUDO rm -f $PREFIX_ROKKO/elemental
  $SUDO ln -s $PREFIX_ROKKO/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION $PREFIX_ROKKO/elemental
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/elementalvars.sh
  $SUDO ln -s $PREFIX_ROKKO/elemental/$build_type/elementalvars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/elementalvars.sh
done

$SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/elementalvars.sh
$SUDO ln -s $PREFIX_ROKKO/elemental/elementalvars.sh $PREFIX_ROKKO/rokkoenv.d/elementalvars.sh
