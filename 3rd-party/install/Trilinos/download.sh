#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

FILE="trilinos-release-$TRILINOS_VERSION_DASHED.tar.gz"
if [ -f $SOURCE_DIR/$FILE ]; then :; else
  wget -O $SOURCE_DIR/$FILE https://github.com/trilinos/Trilinos/archive/$FILE
fi
