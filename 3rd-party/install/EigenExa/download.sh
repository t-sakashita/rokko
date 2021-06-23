#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

FILE="EigenExa-$EIGENEXA_VERSION.tar.gz"
if [ -f $SOURCE_DIR/$FILE ]; then :; else
  wget -O $SOURCE_DIR/$FILE https://www.r-ccs.riken.jp/labs/lpnctrt/projects/eigenexa/$FILE
fi
