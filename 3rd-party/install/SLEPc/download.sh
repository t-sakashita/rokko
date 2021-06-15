#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

FILE="slepc-$SLEPC_VERSION.tar.gz"
if [ -f $SOURCE_DIR/$FILE ]; then :; else
  wget -O $SOURCE_DIR/$FILE http://slepc.upv.es/download/distrib/slepc-$SLEPC_VERSION.tar.gz
fi
