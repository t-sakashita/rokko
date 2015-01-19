#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  $SUDO rm -f $PREFIX_ROKKO/eigenexa
  $SUDO ln -s $PREFIX_ROKKO/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_PATCH_VERSION $PREFIX_ROKKO/eigenexa
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/eigenexa.sh
  $SUDO ln -s $PREFIX_ROKKO/eigenexa/$build_type/eigenexavars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/eigenexavars.sh
done

$SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/eigenexa.sh
$SUDO ln -s $PREFIX_ROKKO/eigenexa/eigenexavars.sh $PREFIX_ROKKO/rokkoenv.d/eigenexavars.sh
