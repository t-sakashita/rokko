#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

cd $BUILD_DIR
rm -rf scalapack-$SCALAPACK_VERSION scalapack-$SCALAPACK_VERSION-build
if test -f $SOURCE_DIR/scalapack-$SCALAPACK_VERSION.tgz; then
  check tar zxf $SOURCE_DIR/scalapack-$SCALAPACK_VERSION.tgz
else
  check wget -O - http://www.netlib.org/scalapack/scalapack-$SCALAPACK_VERSION.tgz | tar zxf -
fi
