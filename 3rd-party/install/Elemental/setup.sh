#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/download.sh

cd $BUILD_DIR
if [ -d Elemental-$ELEMENTAL_VERSION ]; then :; else
  check mkdir -p Elemental-$ELEMENTAL_VERSION
  tar zxf $SOURCE_DIR/Elemental-$ELEMENTAL_VERSION.tar.gz -C Elemental-$ELEMENTAL_VERSION --strip-components=1
  cd Elemental-$ELEMENTAL_VERSION
  if [ -f $SCRIPT_DIR/Elemental-$ELEMENTAL_VERSION.patch ]; then
    echo "applying Elemental-$ELEMENTAL_VERSION.patch"
    check patch -p1 < $SCRIPT_DIR/Elemental-$ELEMENTAL_VERSION.patch
  fi
fi
