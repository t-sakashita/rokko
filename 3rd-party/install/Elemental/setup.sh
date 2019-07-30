#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

cd $BUILD_DIR
rm -rf Elemental-$ELEMENTAL_VERSION*
if [ -f $SOURCE_DIR/Elemental-$ELEMENTAL_VERSION.tgz ]; then
  check tar zxf $SOURCE_DIR/Elemental-$ELEMENTAL_VERSION.tar.gz
else
  check wget -O - https://github.com/elemental/Elemental/archive/v$ELEMENTAL_VERSION.tar.gz | tar zxf -
fi
cd $BUILD_DIR/Elemental-$ELEMENTAL_VERSION
# patch -p1 < $SCRIPT_DIR/elemental-$ELEMENTAL_VERSION.patch
