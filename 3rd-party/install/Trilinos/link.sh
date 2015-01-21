#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  $SUDO rm -f $PREFIX_ROKKO/trilinos
  $SUDO ln -s $PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_PATCH_VERSION $PREFIX_ROKKO/trilinos
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/trilinosvars.sh
  $SUDO ln -s $PREFIX_ROKKO/trilinos/$build_type/trilinosvars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/trilinosvars.sh
done

$SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/trilinosvars.sh
$SUDO ln -s $PREFIX_ROKKO/trilinos/trilinosvars.sh $PREFIX_ROKKO/rokkoenv.d/trilinosvars.sh
