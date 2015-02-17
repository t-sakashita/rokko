#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

cd $BUILD_DIR
rm -rf slepc-$SLEPC_VERSION*
if [ -f $HOME/source/slepc-$SLEPC_VERSION.tar.gz ]; then
  check tar zxf $HOME/source/slepc-$SLEPC_VERSION.tar.gz
else
  check wget -O slepc-$SLEPC_VERSION.tar.gz http://www.grycap.upv.es/slepc/download/download.php?filename=slepc-$SLEPC_VERSION.tar.gz
  check tar zxf slepc-$SLEPC_VERSION.tar.gz
fi
