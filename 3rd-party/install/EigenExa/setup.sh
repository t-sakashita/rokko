#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/download.sh

cd $BUILD_DIR
if [ -d EigenExa-$EIGENEXA_VERSION ]; then :; else
  check mkdir -p EigenExa-$EIGENEXA_VERSION
  tar zxf $SOURCE_DIR/EigenExa-$EIGENEXA_VERSION.tgz -C EigenExa-$EIGENEXA_VERSION --strip-components=1
  cd EigenExa-$EIGENEXA_VERSION
  if [ -f $SCRIPT_DIR/EigenExa-$EIGENEXA_VERSION.patch ]; then
    echo "applying EigenExa-$EIGENEXA_VERSION.patch"
    check patch -p1 < $SCRIPT_DIR/EigenExa-$EIGENEXA_VERSION.patch
  fi
fi
