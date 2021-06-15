#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/download.sh

cd $BUILD_DIR
if [ -f $SOURCE_DIR/trilinos-release-$TRILINOS_VERSION_DASHED.tar.gz ]; then
  check mkdir -p trilinos-$TRILINOS_VERSION
  tar zxf $SOURCE_DIR/trilinos-release-$TRILINOS_VERSION_DASHED.tar.gz -C trilinos-$TRILINOS_VERSION --strip-components=1
  cd trilinos-$TRILINOS_VERSION
  if [ -f $SCRIPT_DIR/trilinos-$TRILINOS_VERSION.patch ]; then
    check patch -p1 < $SCRIPT_DIR/trilinos-$TRILINOS_VERSION.patch
  fi
fi
