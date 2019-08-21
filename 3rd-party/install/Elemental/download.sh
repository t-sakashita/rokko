#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

FILE=Elemental-$ELEMENTAL_VERSION.tar.gz
if [ -f $SOURCE_DIR/$FILE ]; then :; else
  check wget -O $SOURCE_DIR/$FILE https://github.com/elemental/Elemental/archive/v$ELEMENTAL_VERSION.tar.gz
fi
