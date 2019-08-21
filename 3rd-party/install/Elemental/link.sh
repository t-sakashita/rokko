#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  if [ -f $PREFIX_ROKKO/elemental/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION/$build_type/elementalvars.sh ]; then
    $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/elementalvars.sh
    $SUDO ln -s $PREFIX_ROKKO/elemental/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION/$build_type/elementalvars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/elementalvars.sh
  fi
done

if [ -f $PREFIX_ROKKO/elemental/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION/elementalvars.sh ]; then
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/elementalvars.sh
  $SUDO ln -s $PREFIX_ROKKO/elemental/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_RK_REVISION/elementalvars.sh $PREFIX_ROKKO/rokkoenv.d/elementalvars.sh
fi
