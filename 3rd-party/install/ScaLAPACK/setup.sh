#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/download.sh

cd $BUILD_DIR
if [ -d scalapack-$SCALAPACK_VERSION ]; then :; else
  check mkdir -p scalapack-$SCALAPACK_VERSION
  tar zxf $SOURCE_DIR/scalapack-$SCALAPACK_VERSION.tgz -C scalapack-$SCALAPACK_VERSION --strip-components=1
  cd scalapack-$SCALAPACK_VERSION
  if [ -f $SCRIPT_DIR/scalapack-$SCALAPACK_VERSION.patch ]; then
    echo "applying scalapack-$SCALAPACK_VERSION.patch"
    check patch -p1 < $SCRIPT_DIR/scalapack-$SCALAPACK_VERSION.patch
  fi
fi
