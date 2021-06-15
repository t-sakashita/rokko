#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/download.sh

cd $BUILD_DIR
if [ -f $SOURCE_DIR/slepc-$SLEPC_VERSION.tar.gz ]; then
  check mkdir -p slepc-$SLEPC_VERSION
  tar zxf $SOURCE_DIR/slepc-$SLEPC_VERSION.tar.gz -C slepc-$SLEPC_VERSION --strip-components=1
  cd slepc-$SLEPC_VERSION
  if [ -f $SCRIPT_DIR/slepc-$SLEPC_VERSION.patch ]; then
    echo "applying slepc-$SLEPC_VERSION.patch"
    check patch -p1 < $SCRIPT_DIR/slepc-$SLEPC_VERSION.patch
  fi
fi
