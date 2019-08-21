#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/download.sh

cd $BUILD_DIR
if [ -d elpa-$ELPA_VERSION ]; then :; else
  check mkdir -p elpa-$ELPA_VERSION
  tar zxf $SOURCE_DIR/elpa-$ELPA_VERSION.tar.gz -C elpa-$ELPA_VERSION --strip-components=1
  cd elpa-$ELPA_VERSION
  if [ -f $SCRIPT_DIR/elpa-$ELPA_VERSION.patch ]; then
    echo "applying elpa-$ELPA_VERSION.patch"
    check patch -p1 < $SCRIPT_DIR/elpa-$ELPA_VERSION.patch
  fi
fi
