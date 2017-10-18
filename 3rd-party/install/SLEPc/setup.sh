#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

cd $BUILD_DIR
rm -rf slepc-$SLEPC_VERSION*
if [ -f $SOURCE_DIR/slepc-$SLEPC_VERSION.tar.gz ]; then
  check tar zxf $SOURCE_DIR/slepc-$SLEPC_VERSION.tar.gz
else
  check wget http://slepc.upv.es/download/distrib/slepc-$SLEPC_VERSION.tar.gz
  check tar zxf slepc-$SLEPC_VERSION.tar.gz
fi
