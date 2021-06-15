#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

FILE="scalapack-$SCALAPACK_VERSION.tgz"
if [ -f $SOURCE_DIR/$FILE ]; then :; else
  wget -O $SOURCE_DIR/$FILE http://www.netlib.org/scalapack/scalapack-$SCALAPACK_VERSION.tgz
fi
