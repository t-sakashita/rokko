#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_build_dir

cd $BUILD_DIR
rm -rf Elemental-$ELEMENTAL_VERSION*
if [ -f $HOME/source/Elemental-$ELEMENTAL_VERSION.tgz ]; then
  check tar zxf $HOME/source/Elemental-$ELEMENTAL_VERSION.tgz
else
  check wget -O - http://elemental.googlecode.com/files/Elemental-$ELEMENTAL_VERSION.tgz | tar zxf -
fi
cd $BUILD_DIR/Elemental-$ELEMENTAL_VERSION
patch -p1 < $SCRIPT_DIR/elemental-$ELEMENTAL_VERSION.patch
