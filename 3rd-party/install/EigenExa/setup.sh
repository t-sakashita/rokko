#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/download.sh

cd $BUILD_DIR
if [ -d eigenexa-$EIGENEXA_VERSION ]; then :; else
  check mkdir -p eigenexa-$EIGENEXA_VERSION
  tar zxf $SOURCE_DIR/EigenExa-$EIGENEXA_VERSION.tar.gz -C eigenexa-$EIGENEXA_VERSION --strip-components=1
  cd eigenexa-$EIGENEXA_VERSION
  if [ -f $SCRIPT_DIR/EigenExa-$EIGENEXA_VERSION.patch ]; then
    echo "applying EigenExa-$EIGENEXA_VERSION.patch"
    check patch -p1 < $SCRIPT_DIR/EigenExa-$EIGENEXA_VERSION.patch
  fi
fi
