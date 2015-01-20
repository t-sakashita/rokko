#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  $SUDO rm -f $PREFIX_ROKKO/elpa
  $SUDO ln -s $PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_PATCH_VERSION $PREFIX_ROKKO/elpa
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/elpavars.sh
  $SUDO ln -s $PREFIX_ROKKO/elpa/$build_type/elpavars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/elpavars.sh
done

$SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/elpavars.sh
$SUDO ln -s $PREFIX_ROKKO/elpa/elpavars.sh $PREFIX_ROKKO/rokkoenv.d/elpavars.sh
