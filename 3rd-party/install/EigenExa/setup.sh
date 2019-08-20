#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

cd $BUILD_DIR
rm -rf EigenExa-$EIGENEXA_VERSION*
if test -f $SOURCE_DIR/EigenExa-$EIGENEXA_VERSION.tgz; then
  tar zxf $SOURCE_DIR/EigenExa-$EIGENEXA_VERSION.tgz
else
  wget -O $SOURCE_DIR/EigenExa-$EIGENEXA_VERSION.tgz http://www.r-ccs.riken.jp/labs/lpnctrt/assets/img/EigenExa-$EIGENEXA_VERSION.tgz
  tar zxf $SOURCE_DIR/EigenExa-$EIGENEXA_VERSION.tgz
fi
cd $BUILD_DIR/EigenExa-$EIGENEXA_VERSION
patch -p1 < $SCRIPT_DIR/EigenExa-$EIGENEXA_VERSION.patch

