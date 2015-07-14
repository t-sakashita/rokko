#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  $SUDO rm -f $PREFIX_ROKKO/rokko
  $SUDO ln -s $PREFIX_ROKKO/rokko-$ROKKO_VERSION-$ROKKO_RK_REVISION $PREFIX_ROKKO/rokko
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/rokkovars.sh
  $SUDO ln -s $PREFIX_ROKKO/rokko/$build_type/rokkovars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/rokkovars.sh
done

$SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/rokkovars.sh
$SUDO ln -s $PREFIX_ROKKO/rokko/rokkovars.sh $PREFIX_ROKKO/rokkoenv.d/rokkovars.sh
